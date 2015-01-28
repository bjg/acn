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
** gmf: Generate Master Flat file. Multiple Bias files are used to obtain the MEDIAN to produce a master flat
*  all files have each datapoint sorted across each image and the median value is chosen. The program takes a
*  directory as input and assumes all fits files in that directory are flat files to be processed.
*  Flat files can be 2D or 3D. The masterflat output file must be 2D.

*  Paul Doyle 2011, Dublin Institute of Technology V1.0
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
    fprintf(stderr, "Usage: gmf directory outimage \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  gmf ./flatsdir masterflat.fits \n");
}

int main(int argc, char *argv[])
{
    fitsfile **afptr, *outfptr;  /* FITS file pointers */

    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int imagecount = 0,imagecount1=0,counter=0,anaxis, bnaxis, ii;
    long npixels = 1, firstpix[3] = {1,1,1};
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1},cnaxes[3]={1,1,1};
    double *apix, *bpix, *cpix;
    double **dpix; // An array of double pointers to hold the values for pixels across a row and image plane
	long *medianpix;

    // variables to help read list of files
    int count,i,x;
    struct direct **files;
    int file_select();
    int path_max = pathconf(".", _PC_NAME_MAX);
    char fullfilename[path_max];  //to store path and filename
    struct rlimit rl;

	// variables to help with finding the median
	double medianval=0;



    // Verify we have the correct number of parameters

    if (argc != 3) {
        usage();
        exit(0);
    }

    if (argv[1] == NULL) {
        usage();
        bail("Error getting path\n");
    }

    count =  scandir(argv[1], &files, file_select, alphasort);

    // If no files found end the program
    if (count <= 0) {
        bail("No files in this directory\n");
    }

    getrlimit(RLIMIT_NOFILE, &rl);
    if (count > rl.rlim_cur - 3) {
        count = rl.rlim_cur - 3; // stdin, stdout, stderr already open
        fprintf(stderr, "WARNING: Only processing first %d files\n", count);
    }

    afptr = (fitsfile **)calloc(count, sizeof(fitsfile *));

    // Open all of the files
    for (i=0;i<count;++i) {
        snprintf(fullfilename, path_max - 1, "%s%s", argv[1], files[i]->d_name);
        fits_open_file(&afptr[i], fullfilename, READONLY, &status); // open input images
        if (status) {
           fits_report_error(stderr, status); // print error message
           bail(NULL);
        }
    }

    // Use the first file to establish the dimensions for files. All files must have
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
            bail("Error: input images don't have same size\n");
        }

        //calculate the number of images to process
        imagecount1 += bnaxes[2];

    }


	// Allocte space for the median array. The 2D array can hold all the values for each pixel across each
	// image for each column and for a single row. This is done to limit the memory used. We process each
	// row of data across all images first then reuse the storage for each row processed.
	dpix = (double **) malloc(anaxes[1] * sizeof(double *));
    if (dpix == NULL) {
	        bail("Memory allocation error\n");
	}

	if(NULL == dpix){free(dpix); printf("Memory allocation failed while allocating for dpix[].\n"); exit(-1);}

    // Build an array large enough to hold all pixels in each image for a full row
    for (x = 0; x < anaxes[1]; x++) {
    	dpix[x] = (double *) malloc(imagecount1 * sizeof(double)); // memory to hold a single pixel across each image
    	if (dpix[x] == NULL) {
        	bail("Memory allocation error\n");
        }
    }


    // create the new empty output file in the current directory
    if (!fits_create_file(&outfptr, argv[2], &status) ) {

        // Set the image size the same as the images being processed
        cnaxes[0] = anaxes[0];
        cnaxes[1] = anaxes[1];

        fits_create_img(outfptr, LONG_IMG, 2, cnaxes, &status);
        if (status) {
            fits_report_error(stderr, status); // print error message
            bail(NULL);
        }

        npixels = anaxes[0];  // no. of pixels to read in each row

        apix = (double *) malloc(npixels * sizeof(double)); // mem for 1 row to write
        bpix = (double *) malloc(npixels * sizeof(double)); // memory to read each file
        cpix = (double *) malloc(npixels * sizeof(double));  // memory to read each file
        medianpix = (long *) malloc(npixels * sizeof(long));  // memory to read each file

        if (apix == NULL || bpix == NULL || cpix == NULL || medianpix == NULL) {
            bail("Memory allocation error\n");
        }
    } else {
        bail("Output file already exists %s\n",argv[2]);
    }

    // loop over all rows of all of the images
    for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
        // printf ("\nRow %ld",firstpix[1]);
        // Inititalise the apix rows to be zero to start
        bzero((void *) apix, npixels * sizeof(apix[0]));

        if (imagecount > 0) {
            imagecount = 0;
        }

        // Loop through all of the files in our file pointer array
        for (i=0;i<count;++i) {
            fits_get_img_size(afptr[i], 3, bnaxes, &status);  // get the 3D dimension of this file

            // loop over all planes of the file (2D images have 1 plane, 3D images have >1 plane)
            for (firstpix[2] = 1; firstpix[2] <= bnaxes[2]; firstpix[2]++) {
                // Read pixels from images as doubles, regardless of actual datatype.
                // Give starting pixel coordinate and no. of pixels to read.
                // This version does not support undefined pixels in the image.
                if (fits_read_pix(afptr[i], TDOUBLE, firstpix, npixels, NULL, bpix, NULL, &status)) {
                    break;   // jump out of loop on error
                }

                // Add all of the values from each plane/file to the stored array of pixels
                for(ii=0; ii< npixels; ii++) {
                 	dpix[ii][imagecount] = bpix[ii];  // this creates an array of pixels across each image
               }
                imagecount += 1;

            }
        }

        // get the average value for each pixel by dividing the totalled value for all pixels in each
        // data point by the numebr of data points
        //for(ii=0; ii< npixels; ii++) {
        //    cpix[ii] = (apix[ii]/imagecount);
        //}

        for(ii=0; ii< npixels; ii++) {

		     // Sort the pixels into the correct order
		     qsort(dpix[ii],imagecount,sizeof(double),compare_doubles);
		     // Caclulate the MEDIAN in the array
		     // for even numbers of datapoints pick the middle 2 and average them
		     // for odd numebr of datapoints pick the middle element
		     // need to adjust the postion to allow for arrays starting at 0 and not 1
             if ( imagecount % 2 == 0 )
                medianpix[ii]  =(long) (dpix[ii][(imagecount/2)-1] + dpix[ii][(imagecount/2)])/2; // Even
		     else
			    medianpix[ii]  =(long) (dpix[ii][((imagecount+1)/2) -1]); // Odd
		}

		// Write the pixel values out to the output file
		fits_write_pix(outfptr, TLONG, firstpix, npixels, medianpix, &status);

    }  // Move to the next row

		printf("closing files");

    // Close all of the files

    fits_close_file(outfptr, &status);
    for (i=0;i<count;++i) {
        fits_close_file(afptr[i], &status); // close input images
        if (status) {
           fits_report_error(stderr, status); // print error message
           bail(NULL);
        }
    }

    free(afptr);
    free(apix);
    free(bpix);
    free(cpix);
    free(dpix);

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




int compare_doubles (const void *X, const void *Y)
{
       double x = *((double *)X);
       double y = *((double *)Y);

       if (x > y)
       {
               return 1;
       }
       else
       {
               if (x < y)
               {
                       return -1;
               }
               else
               {
                       return 0;
               }
       }
}


int intcmp(const void *v1, const void *v2)
{
  return (*(int *)v1 - *(int *)v2);
}