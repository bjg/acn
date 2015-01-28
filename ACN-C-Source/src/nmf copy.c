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


extern  int alphasort();
int     pr_update_naxis3 ( fitsfile *fptr, int newaxis, int *status);
int 	intcmp(const void *v1, const void *v2);
int 	compare_doubles (const void *X, const void *Y);

/*
** nmf: Normalise Master Flat file. The Master flat file is normalised by dividing each value by the average value of the
*  the data within each pixel. The masterflat output file must be 2D.
*  Paul Doyle 2010, Dublin Institute of Technology
*
*  Modification to allow a flat file to be normalized within a specific subrect. Select an example of an object file
* to create the correct dimensions 
*/

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

void usage(void)
{
    fprintf(stderr, "Usage: nmf biasreducedmasterflat.fits normalisedbiasreducedmasterflat.fits objectfile\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  nmf ./masterflat.fits normalisedmasterflat.fits \n");
}

int main(int argc, char *argv[])
{
    fitsfile *mffptr, *outfptr;  /* FITS file pointers */

    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int imagecount = 0,imagecount1=0,counter=0,anaxis, bnaxis, ii;
    long npixels = 1, firstpix[3] = {1,1,1};
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1},cnaxes[3]={1,1,1};
    long valuecount=0;
    double *apix, sumvalues=0.0,normfactor=0.0;

    // variables to help read list of files
    int count,i,x;

	// variables to help with finding the median
	double medianval=0;

    // Verify we have the correct number of parameters

    if (argc != 3) {
        usage();
        exit(0);
    }

    if (argv[1] == NULL) {
        usage();
        bail("Error getting input file \n");
    }

    if (argv[2] == NULL) {
	        usage();
	        bail("Error getting output file \n");
    }

    // Open the master flat file
    fits_open_file(&mffptr, argv[1], READONLY, &status); // open input images
    if (status) {
       fits_report_error(stderr, status); // print error message
       bail(NULL);
    }


    // Read the input file to establish the dimensions for files.

    fits_get_img_dim(mffptr, &anaxis, &status);  // read dimensions of the file
    fits_get_img_size(mffptr, 3, anaxes, &status);

    if (status) {
        fits_report_error(stderr, status); // print error message
        bail(NULL);
    }

    if (anaxis > 2) {
       bail("Error: Master Flat File %s in an images with > 2 dimensions and is not supported\n",argv[1]);
    }


    // create the new empty output file in the current directory
    if (!fits_create_file(&outfptr, argv[2], &status) ) {

        // Set the image size the same as the images being processed
        cnaxes[0] = anaxes[0];
        cnaxes[1] = anaxes[1];

        fits_create_img(outfptr, DOUBLE_IMG, 2, cnaxes, &status);
        if (status) {
            fits_report_error(stderr, status); // print error message
            bail(NULL);
        }

        npixels = anaxes[0];  // no. of pixels to read in each row

        apix = (double *) malloc(npixels * sizeof(double)); // mem for 1 row to write

        if (apix == NULL) {
            bail("Memory allocation error\n");
        }
    } else {
        bail("Output file already exists %s\n",argv[2]);
    }


    // loop over all rows of all of the images
    // calculate the average value for pixels
    for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
        // Inititalise the apix rows to be zero to start
        bzero((void *) apix, npixels * sizeof(apix[0]));

        if (fits_read_pix(mffptr, TDOUBLE, firstpix, npixels, NULL, apix, NULL, &status)) {
              printf("Failed to read row %ld \n",firstpix[1]);
              break;   // jump out of loop on error
		}

        // get the average value for each pixel by dividing the totalled value for all pixels in each
        // data point by the numebr of data points
    	for(ii=0; ii< npixels; ii++) {
          	sumvalues += apix[ii];
          	valuecount++;
    	}
	}

	normfactor = sumvalues/valuecount;
	printf ("finished reading %ld values total = %lf average = %lf\n",valuecount,sumvalues,normfactor);

    // loop over all rows again subtracting the bias values and dividing by the normalisation factor normfactor
    // calculate the average value for pixels
    for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
        // Inititalise the apix rows to be zero to start
        bzero((void *) apix, npixels * sizeof(apix[0]));

        if (fits_read_pix(mffptr, TDOUBLE, firstpix, npixels, NULL, apix, NULL, &status)) {
              bail("Failed to read Flat File row %ld \n",firstpix[1]);
		}

        // get the average value for each pixel by dividing the totalled value for all pixels in each
        // data point by the numebr of data points
    	for(ii=0; ii< npixels; ii++) {
			apix[ii] = apix[ii]/normfactor;
    	}

        // write the normalised values to the output file
    	fits_write_pix(outfptr, TDOUBLE, firstpix, npixels, apix, &status);

	}

    // Close all of the files

    fits_close_file(outfptr, &status);
    fits_close_file(mffptr, &status);
	if (status) {
           fits_report_error(stderr, status); // print error message
           bail(NULL);
    }

    free(apix);

    exit(0);
}
