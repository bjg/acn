#include <string.h>
#include "fitsio.h"
#include <sys/types.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/resource.h>
#include <strings.h>

// These are used to determine which message should be printed when debugging
#define  DEBUGLEVEL1 1
#define  DEBUGLEVEL2 2
#define  DEBUGLEVEL3 3

extern  int alphasort();
int     pr_update_naxis3 ( fitsfile *fptr, int newaxis, int *status);

/*
** gmb: Generate Master Bias file. Multiple Bias files are combined with their data values
*  averaged across all images to produce a master bias file.
*  All files have each datapoint added and divided by the number of values. The program takes a
*  directory as input and assumes all fits files in that directory are bias files to be processed.
*  Bias files can be 2D or 3D. The masterbias output file must be 2D.
*
*  Author/Copyright: Paul Doyle 2011 V1.0
*/

//
// Terminate the program with a specified effort message
//
void bail(const char *msg, ...)
{
    va_list arg_ptr;

    va_start(arg_ptr, msg);
    if (msg) {
        vfprintf(stderr, msg, arg_ptr);
    }
    va_end(arg_ptr);
    fprintf(stderr, "\nAborting...\n");

    exit(1);
}

// Allow debug messages to be sent to the screen
// It will only print out debug messages that are at
// of the same or greater than the requested debug level
// If a user debug level = 0 the don't print debug messages
// which are Level 1 or greather. The higher the debug message
// the more deailed the output is likely to be

void debug(int userlevel, int debuglevel, const char *msg, ...)
{
	if (userlevel >= debuglevel > 0) {
		va_list arg_ptr;

    	va_start(arg_ptr, msg);
    	if (msg) {
    	    vfprintf(stderr, msg, arg_ptr);
    	}
    	va_end(arg_ptr);
	}
}

//
// Inform the user of the correct use of this program
//
void usage(void)
{
    fprintf(stderr, "This program will generate a Master Bias FITS file from\n");
    fprintf(stderr, "a directory of bias files\n\n");
    fprintf(stderr, "You can optionally run this program in debug mode for extra output\n");
    fprintf(stderr, "Usage: gmb directory outimage {-debug1|-debug2|-debug3} \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  gmb ./biasdir masterbias.fits \n");
    fprintf(stderr, "  gmb ./biasdir masterbias.fits -debug1\n");
}

int main(int argc, char *argv[])
{

    // Variable to help process the fits files
    fitsfile **afptr, *outfptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int imagecount = 0, anaxis, bnaxis, ii;
    int debuglevel =0; /* Starting debug level is off */
    long npixels = 1, firstpix[3] = {1,1,1};
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1},cnaxes[3]={1,1,1};
    double *apix, *bpix;
    double *cpix;

    // variables to help read list of files
    int count,i;
    struct direct **files;
    int file_select();
    int path_max = pathconf(".", _PC_NAME_MAX);
    char fullfilename[path_max];  //to store path and filename
    struct rlimit rl;

    // Verify we have the correct number of parameters
    // If there are 3 parameters, then the only valid ones are
    // d1, d2, d3. Anything else causes program to stop
    if (argc    == 4)    {
		if (strcmp(argv[3],"-debug1")== 0 ) debuglevel=1;
		if (strcmp(argv[3],"-debug2")== 0)  debuglevel=2;
		if (strcmp(argv[3],"-debug3") ==0)  debuglevel=3;
		if (debuglevel==0) {usage(); exit(0);}

	} else {
		if (argc    != 3)    {  usage(); exit(0);}
	}

    // count the number of files in the directory to process
    // If no files found end the program
    count =  scandir(argv[1], &files, file_select, alphasort);
    if (count <= 0) {  bail("No files to process in this directory\n");  }

	// If we cannot open all of the files, inform the user
	// The limit to the number of files openable is set in the OS
	getrlimit(RLIMIT_NOFILE, &rl);
    if (count > rl.rlim_cur - 3) {
        count = rl.rlim_cur - 3; // stdin, stdout, stderr already open
        fprintf(stderr, "WARNING: Only processing first %d files\n, increase RLIMIT", count);
    }

    // Allocate enough memory for an array of filepointers
 	afptr = (fitsfile **)calloc(count, sizeof(fitsfile *));
    debug(debuglevel,DEBUGLEVEL1,"Number of .fits files = %d\n",count);

    // Open all of the files identified
    for (i=0;i<count;++i) {

		// Generate the complete file name to open
        snprintf(fullfilename, path_max - 1, "%s%s", argv[1], files[i]->d_name);
        fits_open_file(&afptr[i], fullfilename, READONLY, &status); // open input images
        if (status) {
           fits_report_error(stderr, status); // print error message
           bail("failed to open an input file");
        }
    }

    // Use the first file to establish the dimensions for images. All files must have
    // same basic size for a  image, but there may have multiple images in a file.
    fits_get_img_dim(afptr[0], &anaxis, &status);  // read dimensions of each file
    fits_get_img_size(afptr[0], 3, anaxes, &status);
    if (status) {
        fits_report_error(stderr, status); // print error message
        bail(NULL);
    }

    if (anaxis > 3) {
       bail("Error: File %s in an images with > 3 dimensions and is not supported\n",files[0]->d_name);
    }

    // additional files must have the same reference size for an image as the first file for processing
    // to make sense.
    for (i=0;i<count;++i) {
        fits_get_img_dim(afptr[i], &bnaxis, &status);  // read dimensions of each file
        fits_get_img_size(afptr[i], 3, bnaxes, &status);

        if (status) {
            fits_report_error(stderr, status); // print error message
            bail(NULL);
        }

        if (bnaxis > 3) {
            bail("Error: File %s in an images with > 3 dimensions and is not supported\n",files[i]->d_name);
        }

        // We only need to check the image size, not the number of images in a file.
        if (( anaxes[0] != bnaxes[0] || anaxes[1] != bnaxes[1] )) {
            bail("Error: File %s input image is of a different size don't have same size\n",files[i]->d_name);
        }
    }


    // create the new empty output file in the current directory
    // Even if the bias files being read contained stacked images
    // the output file must be a 2D image only.
    // File must not already exist and it will not overwrite a file that does
    // already exist.
    if (!fits_create_file(&outfptr, argv[2], &status) ) {

        // Set the image size the same as the images being processed
        cnaxes[0] = anaxes[0];
        cnaxes[1] = anaxes[1];

        fits_create_img(outfptr, DOUBLE_IMG, 2, cnaxes, &status);
        if (status) {
            fits_report_error(stderr, status); // print error message
            bail(NULL);
        }

		// The memory requirements is based processing one row at a time.
		// These arrays allow us to store and calculate the average value of each pixel
		// across all images.
		//
		// apix contains the sum of the values for each image
		// bpix contains the values read for an image
		// cpix contains the average value for the pixels across all images
        npixels = anaxes[0];  // no. of pixels to read in each row
        apix = (double *) malloc(npixels * sizeof(double)); // mem for 1 row to write
        bpix = (double *) malloc(npixels * sizeof(double)); // memory to read each file
        cpix = (double *) malloc(npixels * sizeof(double)); // memory to read each file

        if (apix == NULL || bpix == NULL || cpix == NULL) {
            bail("Memory allocation error\n");
        }
    } else {
        bail("Output file already exists %s\n",argv[2]);
    }



	// This is the main processing loop.
	// The process is to loop over each row from the top to the bottom of the image
	// For each row open all the files, and within each file loop through all the
	// image if the file contains stacked images.
	//


	// Loop through each row, from top of image to bottom of image
    for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
    	debug(debuglevel,DEBUGLEVEL1,"Processing Row = %ld\n",firstpix[1]);

        // Inititalise the apix rows to be zero to start
        bzero((void *) apix, npixels * sizeof(apix[0]));

        // We need to know the number of images so we can correctly
        // calculate the average value for each pixel. We rest the
        // image count to 0 when we start a new row.
        if (imagecount > 0) {
            imagecount = 0;
        }

        // Loop through all of the files in our file pointer array, and in each
        // file loop through the image.
        for (i=0;i<count;++i) {
            fits_get_img_size(afptr[i], 3, bnaxes, &status);  // get the dimension of this file
    		debug(debuglevel,DEBUGLEVEL2,"Processing File = %s\n",files[i]->d_name);

            // loop over all planes of the file (2D images have 1 plane, 3D images have >1 plane)
            for (firstpix[2] = 1; firstpix[2] <= bnaxes[2]; firstpix[2]++) {
                // Read pixels from images as doubles, regardless of actual datatype.
                // Give starting pixel coordinate and no. of pixels to read.
                // This version does not support undefined pixels in the image.
                if (fits_read_pix(afptr[i], TDOUBLE, firstpix, npixels, NULL, bpix, NULL, &status)) {
                    break;   // jump out of loop on error
                }

    			debug(debuglevel,DEBUGLEVEL3,"Processing Image = %ld\n",imagecount+1);

                // add the values from the current image row to our apix array. As we do
                // this for all images we obtain the total of all values for the pixels
                // within the row.
                for(ii=0; ii< npixels; ii++) {
                    apix[ii] += bpix[ii];
                }
                imagecount += 1;
            }
        }

        // For the current row, all images within all files have been added to apix
        // now we get the average value for each pixel by dividing the totalled value for all pixels in each
        // data point by the numebr of data points
        for(ii=0; ii< npixels; ii++) {
            cpix[ii] = (apix[ii]/imagecount);
        }

        // Write the pixel values out to the output file
        fits_write_pix(outfptr, TDOUBLE, firstpix, npixels, cpix, &status);

    }  // Move to the next row


    // Close all of the files
    fits_close_file(outfptr, &status);
    for (i=0;i<count;++i) {
        fits_close_file(afptr[i], &status); // open input images
        if (status) {
           fits_report_error(stderr, status); // print error message
           bail(NULL);
        }
    }

	// Free all of the memory allocated
    free(afptr);
    free(apix);
    free(bpix);
    free(cpix);

    exit(0);
}

int file_select(struct direct   *entry) {

    char *ptr;

    if ((strcmp(entry->d_name, ".")== 0) ||    (strcmp(entry->d_name, "..") == 0))
        return (FALSE);

    /* Check for filename extensions */
    ptr = strrchr(entry->d_name, '.');
    return ((ptr != NULL) && (strcmp(ptr, ".fits") == 0));
}

