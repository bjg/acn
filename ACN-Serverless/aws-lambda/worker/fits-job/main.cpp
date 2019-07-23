// main.cpp
#include <aws/core/Aws.h>
#include <aws/core/utils/logging/LogLevel.h>
#include <aws/core/utils/logging/ConsoleLogSystem.h>
#include <aws/core/utils/logging/LogMacros.h>
#include <aws/core/utils/json/JsonSerializer.h>
#include <aws/core/utils/HashingUtils.h>
#include <aws/core/platform/Environment.h>
#include <aws/core/client/ClientConfiguration.h>
#include <aws/core/auth/AWSCredentialsProvider.h>
#include <aws/s3/S3Client.h>
#include <aws/s3/model/GetObjectRequest.h>
#include <aws/lambda-runtime/runtime.h>
#include <iostream>
#include <memory>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>  /* for open flags */
#include <limits.h> /* for PATH_MAX */
#include <unistd.h>
#include <fitsio.h>

class TempFile
{
private:
    char *fname;
    int fd;
    size_t len;

public:
    TempFile()
    {
        static char const *tmpl = "/tmp/fitsfileXXXXXX";
        fname = new char[PATH_MAX];
        strcpy(fname, tmpl); /* Copy template */
        fd = mkstemp(fname); /* Create and open temp file */
        if (fd == -1)
        {
            throw strerror(errno);
        }
        len = 0;
    }

    size_t put(unsigned char *p, size_t n)
    {
        size_t count = write(fd, p, n);
        len += count;
        return count;
    }

    int flush()
    {
        int status = close(fd);
        fd = -1;
        return status;
    }

    const char *name() { return (const char *)fname; }

    const size_t size() { return len; }

    ~TempFile()
    {
        if (fd != -1)
            close(fd); /* Close file (if not already so) */
        unlink(fname); /* Remove it */
        delete fname;
    }
};

using namespace aws::lambda_runtime;

int download_and_decode_fits_file(
    Aws::S3::S3Client const &client,
    Aws::String const &bucket,
    Aws::String const &key,
    Aws::String &decoded_output);

Aws::String decode_header(Aws::String const &filename, Aws::String &output);
char const TAG[] = "LAMBDA_ALLOC";

// Example:
// '{\"s3bucket\":\"tudublincs-research-acn-images\",\"s3key\":\"WFPC2u5780205r_c0fx.fits\" }'

static invocation_response my_handler(invocation_request const &req, Aws::S3::S3Client const &client)
{
    using namespace Aws::Utils::Json;
    JsonValue json(req.payload);
    if (!json.WasParseSuccessful())
    {
        return invocation_response::failure("Failed to parse input JSON", "InvalidJSON");
    }

    auto v = json.View();

    if (!v.ValueExists("s3bucket") || !v.ValueExists("s3key") || !v.GetObject("s3bucket").IsString() ||
        !v.GetObject("s3key").IsString())
    {
        return invocation_response::failure("Missing input value s3bucket or s3key", "InvalidJSON");
    }

    auto bucket = v.GetString("s3bucket");
    auto key = v.GetString("s3key");

    AWS_LOGSTREAM_INFO(TAG, "Attempting to download file from s3://" << bucket << "/" << key);

    Aws::String header_content;
    auto status = download_and_decode_fits_file(client, bucket, key, header_content);
    if (status == -1)
    {
        return invocation_response::failure(header_content, "DownloadFailure");
    }

    return invocation_response::success(header_content, "text/plain");
}

std::function<std::shared_ptr<Aws::Utils::Logging::LogSystemInterface>()> GetConsoleLoggerFactory()
{
    return [] {
        return Aws::MakeShared<Aws::Utils::Logging::ConsoleLogSystem>(
            "console_logger", Aws::Utils::Logging::LogLevel::Trace);
    };
}

// Create a temporary local copy of the download object to satisfy the FITS API requirements
TempFile *download_file(Aws::IOStream &stream, Aws::String &output)
{
    stream.seekg(0, stream.beg);
    auto *tf = new TempFile();

    char streamBuffer[1024 * 4];
    while (stream.good())
    {
        stream.read(streamBuffer, sizeof(streamBuffer));
        auto bytesRead = stream.gcount();

        if (bytesRead > 0)
        {
            tf->put((unsigned char *)streamBuffer, bytesRead);
        }
    }
    tf->flush();

    return tf;
}

Aws::String decode_header(const char *fname)
{
    char card[FLEN_CARD]; /* Standard string lengths defined in fitsio.h */
    int status = 0;
    int single = 0, hdupos, nkeys, ii;
    Aws::String output("");
    fitsfile *fptr;

    if (!fits_open_file(&fptr, fname, READONLY, &status))
    {
        fits_get_hdu_num(fptr, &hdupos); /* Get the current HDU position */

        /* List only a single header if a specific extension was given */
        if (hdupos != 1 || strchr(fname, '['))
            single = 1;

        for (; !status; hdupos++) /* Main loop through each extension */
        {
            fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */

            output.append("Header listing for HDU ").append(std::to_string(hdupos)).append("\n");

            for (ii = 1; ii <= nkeys; ii++)
            { /* Read and print each keywords */

                if (fits_read_record(fptr, ii, card, &status))
                    break;
                output.append(card).append("\n");
            }
            output.append("END").append("\n\n"); /* terminate listing with END */

            if (single)
                break; /* quit if only listing a single header */

            fits_movrel_hdu(fptr, 1, NULL, &status); /* try to move to next HDU */
        }

        if (status == END_OF_FILE)
            status = 0; /* Reset after normal error */

        fits_close_file(fptr, &status);
    }

    return output;
}

int download_and_decode_fits_file(
    Aws::S3::S3Client const &client,
    Aws::String const &bucket,
    Aws::String const &key,
    Aws::String &decoded_output)
{
    using namespace Aws;

    S3::Model::GetObjectRequest request;
    request.WithBucket(bucket).WithKey(key);

    auto outcome = client.GetObject(request);
    if (outcome.IsSuccess())
    {
        Aws::String output("");

        AWS_LOGSTREAM_INFO(TAG, "Download completed!");
        auto &s = outcome.GetResult().GetBody();
        TempFile *fits_file = download_file(s, decoded_output);
        decoded_output = decode_header(fits_file->name());
        delete fits_file;
        return 0;
    }
    else
    {
        AWS_LOGSTREAM_ERROR(TAG, "Failed with error: " << outcome.GetError());
        decoded_output = outcome.GetError().GetMessage();
        return -1;
    }
}

int main()
{
    using namespace Aws;
    SDKOptions options;
    options.loggingOptions.logLevel = Aws::Utils::Logging::LogLevel::Trace;
    options.loggingOptions.logger_create_fn = GetConsoleLoggerFactory();
    InitAPI(options);
    {
        Client::ClientConfiguration config;
        config.region = Aws::Environment::GetEnv("AWS_REGION");
        config.caFile = "/etc/pki/tls/certs/ca-bundle.crt";

        auto credentialsProvider = Aws::MakeShared<Aws::Auth::EnvironmentAWSCredentialsProvider>(TAG);
        S3::S3Client client(credentialsProvider, config);
        auto handler_fn = [&client](aws::lambda_runtime::invocation_request const &req) {
            return my_handler(req, client);
        };
        run_handler(handler_fn);
    }
    ShutdownAPI(options);
    return 0;
}
