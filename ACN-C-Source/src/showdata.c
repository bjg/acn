#include <string.h>
#include <stdio.h>
#include "fitsio.h"

/*
** Print entire file or points within in a 2D images 
*/

int main(int argc, char *argv[])
{
    fitsfile *afptr ;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int atype, btype, anaxis, bnaxis, check = 1, ii, op,errorcounter=0;
    long npixels = 1, firstpix[3] = {1,1,1}, ntodo;
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1};
    double *apix, value;

    if (argc != 3) {
      printf("Usage: showdata image1 plane x y\n");
      printf("\n");
      printf("Examples: show first pixel in first plane\n");
      printf("  showdata in1.fits 1 1 1   \n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input images */
    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_size(afptr, 3, anaxes, &status);

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    if (anaxis > 3) {
       printf("Error: images with > 3 dimensions are not supported\n");
       check = 0;
    }


      npixels = anaxes[0];  /* no. of pixels to read in each row */
      apix = (double *) malloc(npixels * sizeof(double)); /* mem for 1 row */

      if (apix == NULL) {
        printf("Memory allocation error\n");
        return(1);
  	  }


  /* loop over all planes of the cube (2D images have 1 plane) */
      for (firstpix[2] = 1; firstpix[2] <= 2; firstpix[2]++)
      {
        /* loop over all rows of the plane */
        for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++)
        {
          /* Read both images as doubles, regardless of actual datatype.  */
          /* Give starting pixel coordinate and no. of pixels to read.    */
          /* This version does not support undefined pixels in the image. */

          if (fits_read_pix(afptr, TDOUBLE, firstpix, npixels, NULL, apix,
                            NULL, &status))
	    	break;   /* jump out of loop on error */
          for(ii=0; ii< npixels; ii++)
	  	printf ("plane = %ld, r %ld, c %d, %*.*f\n",firstpix[2],firstpix[1],ii,11,10,apix[ii]);
			  	errorcounter++;
        }
      }    /* end of loop over planes */


    printf ("The number of differences = %d\n",errorcounter);

    free(apix);

    fits_close_file(afptr, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}

