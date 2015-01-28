#include <string.h>
#include <stdio.h>
#include "fitsio.h"

/*
** Compare 2 2D images for exactly the same values and print out any differences
*/

int main(int argc, char *argv[])
{
    fitsfile *afptr, *bfptr, *outfptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int atype, btype, anaxis, bnaxis, check = 1, ii, op,counter=0,errorcounter=0;
    long npixels = 1, firstpix[3] = {1,1,1}, ntodo;
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1};
    double *apix, *bpix, value;

    if (argc != 3) {
      printf("Usage: imarith image1 image2 \n");
      printf("\n");
      printf("Examples: \n");
      printf("  compare in1.fits in2.fits  compare the 2 files\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input images */
    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }
    fits_open_file(&bfptr, argv[2], READONLY, &status);
    if (status) {
	       fits_report_error(stderr, status); /* print error message */
			return(status);
    }

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_dim(bfptr, &bnaxis, &status);
    fits_get_img_size(afptr, 3, anaxes, &status);
    fits_get_img_size(bfptr, 3, bnaxes, &status);

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    if (anaxis > 3) {
       printf("Error: images with > 3 dimensions are not supported\n");
       check = 0;
    }
         /* check that the input 2 images have the same size */
    else if ( ( anaxes[0] != bnaxes[0] ||
			  anaxes[1] != bnaxes[1] ||
			  anaxes[2] != bnaxes[2] ) ) {
        printf("Error: input images don't have same size\n");
	printf ("anaxes[1] = %ld, bnaxeses[1] = %ld \n",anaxes[1],bnaxes[1]);
	printf ("anaxes[2] = %ld, bnaxeses[2] = %ld \n",anaxes[2],bnaxes[2]);
       //exit(0);
    }


      npixels = anaxes[0];  /* no. of pixels to read in each row */
      apix = (double *) malloc(npixels * sizeof(double)); /* mem for 1 row */
      bpix = (double *) malloc(npixels * sizeof(double));

      if (apix == NULL || bpix == NULL) {
        printf("Memory allocation error\n");
        return(1);
  	  }


  /* loop over all planes of the cube (2D images have 1 plane) */
      for (firstpix[2] = 1; firstpix[2] <= 1; firstpix[2]++)
     // for (firstpix[2] = 1; firstpix[2] <= anaxes[2]; firstpix[2]++)
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
	  	  if (fits_read_pix(bfptr, TDOUBLE, firstpix, npixels,
				      NULL, bpix, NULL, &status))
	    	break;   /* jump out of loop on error */

          for(ii=0; ii< npixels; ii++) 
	      if ((apix[ii]-bpix[ii] > 0.001) || (bpix[ii]-apix[ii] > 0.001) ) {
			  	printf ("diff = %f   ",apix[ii]-bpix[ii]);
				printf ("plane = %ld, r %ld, c %d, %*.*f, %*.*f\n",firstpix[2],firstpix[1],ii,11,14,apix[ii],11,14,bpix[ii]);
			  	errorcounter++;
	  
		}
		else counter++;
        }
      }    /* end of loop over planes */


    printf ("The number of items checked %d - differences = %d\n",counter,errorcounter);

    free(apix);
    free(bpix);

    fits_close_file(afptr, &status);
    fits_close_file(bfptr, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}

