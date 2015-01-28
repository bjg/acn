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
** bmf: Bias reduce Master Flat file. The Master flat file has the master bias remvoved from each pixel.
*  Master Flat file should be 2D. The masterflat output file must be 2D. The Master Bias and MasterFlat must
*  have the same dimensions
*  Paul Doyle 2010, Dublin Institute of Technology
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
    fprintf(stderr, "Usage: bmf masterflat.fits  masterbias.fits output.fits \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  bmf ./masterflat.fits ./masterbias.fits biasreducemasterflat.fits \n");
}

int main(int argc, char *argv[])
{
    fitsfile *mffptr, *mbfptr, *outfptr;  /* FITS file pointers */

    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int imagecount = 0,imagecount1=0,counter=0,anaxis, bnaxis, ii;
    long npixels = 1, firstpix[3] = {1,1,1};
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1},cnaxes[3]={1,1,1};
    double *apix,*bpix, valuecount=0.0,sumvalues=0.0,normfactor=0.0;

    // variables to help read list of files
    int count,i,x;

	// variables to help with finding the median
	double medianval=0;

    // Verify we have the correct number of parameters

    if (argc != 4) {
        usage();
        exit(0);
    }

    if (argv[1] == NULL) {
        usage();
        bail("Error getting input file \n");
    }

    if (argv[2] == NULL) {
        usage();
        bail("Error getting masterbias file \n");
    }

    if (argv[3] == NULL) {
	        usage();
	        bail("Error getting output file \n");
    }

    // Open the master flat file
    fits_open_file(&mffptr, argv[1], READONLY, &status); // open input images
    if (status) {
       fits_report_error(stderr, status); // print error message
       bail(NULL);
    }

    // Open the master bias file
    fits_open_file(&mbfptr, argv[2], READONLY, &status); // open input images
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

	// Verify that the Master Bias and Master Flat are the same dimensions
    fits_get_img_dim(mbfptr, &bnaxis, &status);  // read dimensions of each file
    fits_get_img_size(mbfptr, 3, bnaxes, &status);

    if (status) {
         fits_report_error(stderr, status); // print error message
         bail(NULL);
    }

    if (bnaxis > 2)
         bail("Error: Master Bias File %s in an images with > 2 dimensions and is not supported\n",argv[2]);

    // We only need to check the image size, not the number of images in a file.
    if (( anaxes[0] != bnaxes[0] || anaxes[1] != bnaxes[1] ))
        bail("Error: input images don't have same size\n");


    // create the new empty output file in the current directory
    if (!fits_create_file(&outfptr, argv[3], &status) ) {

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
        bpix = (double *) malloc(npixels * sizeof(double)); // mem for 1 row to write

        if (apix == NULL) {
            bail("Memory allocation error\n");
        }
    } else {
        bail("Output file already exists %s\n",argv[3]);
    }


    // loop over all rows of all of the images
    // and remove the bias from them
    for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
        // Inititalise the apix rows to be zero to start
        bzero((void *) apix, npixels * sizeof(apix[0]));
        bzero((void *) bpix, npixels * sizeof(bpix[0]));

        if (fits_read_pix(mffptr, TDOUBLE, firstpix, npixels, NULL, apix, NULL, &status)) {
              bail("Failed to read Flat File row %ld \n",firstpix[1]);
		}

        if (fits_read_pix(mbfptr, TDOUBLE, firstpix, npixels, NULL, bpix, NULL, &status)) {
              bail("Failed to read Bias row %ld \n",firstpix[1]);
		}

        // remove the bias values from the master flat
    	for(ii=0; ii< npixels; ii++) {
          	apix[ii] -= bpix[ii];
    	}

    	// write the bias reduced values to the output file
		fits_write_pix(outfptr, TDOUBLE, firstpix, npixels, apix, &status);

	}


    // Close all of the files

    fits_close_file(outfptr, &status);
    fits_close_file(mffptr,  &status);
    fits_close_file(mbfptr,  &status);

	if (status) {
           fits_report_error(stderr, status); // print error message
           bail(NULL);
    }

    free(apix);
    free(bpix);

    exit(0);
}
