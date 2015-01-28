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

double 	pr_julian_date (int year, int month, int day,int hour, int minute, double second);
int 	pr_update_date ( fitsfile *fptr, double jd, int *status);


extern  int alphasort();
int     pr_update_naxis3 ( fitsfile *fptr, int newaxis, int *status);
int 	intcmp(const void *v1, const void *v2);
int 	compare_doubles (const void *X, const void *Y);

/*
** rrf: Reduce Raw object file.
*       The requirement is to read in a series of raw object files from a directory and perform the
*		following operation per pixel
*       (rawdata-bias value)/normalisedMaster Flat value
*       In some cases the raw object file may only use a subrectangle so the data dimensions will be
*       potentiall different.
*		The processed files are then output to a specified directory
*
*		Paul Doyle 2010, Dublin Institute of Technology
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
    fprintf(stderr, "Usage: rrf directory masterflat masterbias outdir \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  cleanobjectfile ./sourcedir ./masterflat.fits ./masterbias.fits ./outdir \n");
}

int main(int argc, char *argv[])
{
    fitsfile *datafptr, *mffptr, *mbfptr, *outfptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
	int anaxis, bnaxis,cnaxis, ii;
	long cval =0;

    int imagecount = 0,imagecount1=0,counter=0,subrectdim[4] = {1,1,1,1};
    long npixels = 1, ndpixels=1, firstpix[3] = {1,1,1},indexpix[3] = {1,1,1};
    long anaxes[3] = {1,1,1}, bnaxes[3]={1,1,1},cnaxes[3]={1,1,1};
    double *apix,*bpix, *cpix, valuecount=0.0,sumvalues=0.0,normfactor=0.0;

    // variables to help read list of files
    int count,i,x,subrect=0;
    struct direct **files;
    int file_select();
    int path_max = pathconf(".", _PC_NAME_MAX);
    char fullfilename[path_max];  //to store path and filename
    struct rlimit rl;

	int length;
	 char *to,*subrectstring;
	// variables to help with finding the median
	double medianval=0;

	char date[30],keyval[FLEN_VALUE],comment[FLEN_COMMENT];
    int	year,month,day,hour,minute;
 	double second,jd;
 	char iskey[21];

 	int ic =0, rc=0;

    // Verify we have the correct number of parameters
    //
    // There are 5 parameters to verify, the first is always the program name
    // A source directory should be used to read files for processing but a
    // destination directory is also required
    // The name of the output file will be "BF"_originalName.fits
	//

    if (argc != 5)       { usage();  exit(0); }
    if (argv[1] == NULL) { usage();bail("Error getting Source directory \n");  }
    if (argv[2] == NULL) { usage();bail("Error getting Master flat file \n");  }
    if (argv[3] == NULL) { usage();bail("Error getting Master bias file \n");  }
    if (argv[4] == NULL) { usage();bail("Error getting Output directory \n");  }

    count =  scandir(argv[1], &files, file_select, alphasort);

    printf ("Processing %d files \n",count);

	// Next we process the Flat and Bias files which must pass a series
	// of checks.
	// 1. We should be able to open the 2 files
	// 2. The dimensions should be 2D for both files
	// 3. The 2D dimensions should be identical for both files.

    // Open the master flat file
    fits_open_file(&mffptr, argv[2], READONLY, &status); // open input images
    if (status) {
       fits_report_error(stderr, status); // print error message
       bail(NULL);
    }

    // Open the master bias file
    fits_open_file(&mbfptr, argv[3], READONLY, &status); // open input images
    if (status) {
       fits_report_error(stderr, status); // print error message
       bail(NULL);
    }

    // Read the flat and bias files to establish the dimensions are identical.
    fits_get_img_dim(mffptr, &anaxis, &status);  // read dimensions of the file
    fits_get_img_size(mffptr, 3, anaxes, &status);
    if (status) {
        fits_report_error(stderr, status); // print error message
        bail(NULL);
    }
    if (anaxis > 2)
       bail("Error: Master Flat File %s in an images with > 2 dimensions and is not supported\n",argv[1]);

	// Verify that the Master Bias and Master Flat are the same dimensions
    fits_get_img_dim(mbfptr, &bnaxis, &status);  // read dimensions of each file
    fits_get_img_size(mbfptr, 3, bnaxes, &status);

    if (status) {
         fits_report_error(stderr, status); // print error message
         bail(NULL);
    }
    if (bnaxis > 2)
         bail("Error: Master Bias File %s in an images with > 2 dimensions and is not supported\n",argv[2]);

    // Bias and Master files should be the same size.
    if (( anaxes[0] != bnaxes[0] || anaxes[1] != bnaxes[1] ))
        bail("Error: input images don't have same size\n");



	//
	// Process each of the data files, Remembering that the data files may be Cubed files
	//
    for (i=0;i<count;++i) {

			// Open the input files
        	snprintf(fullfilename, path_max - 1, "%s%s", argv[1], files[i]->d_name);
        	fits_open_file(&datafptr, fullfilename, READONLY, &status); // open input images
        	if (status) {
        	   fits_report_error(stderr, status); // print error message
        	   bail(NULL);
        	}

			// Check the dimension of the DATA file */
			// cnaxis give the dimensions */
			fits_get_img_dim(datafptr, &cnaxis, &status);  // read dimensions

			// Next we get the dimension filled in our 3D array anaxes
			fits_get_img_size(datafptr, 3, cnaxes, &status);
			if (status) {
			   fits_report_error(stderr, status); /* print error message */
			   return(status);
			}


			// Create the output file with the same dimensions as the data file.
			//
			snprintf(fullfilename, path_max - 1, "%s/BF_%s", argv[4], files[i]->d_name);
			printf ("Writing ...%s\n",fullfilename);

			if (!fits_create_file(&outfptr, fullfilename, &status) ) {
			   // copy all the header keywords from first image to new output file
			   fits_copy_header(datafptr, outfptr, &status);
			   ndpixels = cnaxes[0];  // no. of data pixels to read in each row

			   apix = (double *) malloc(ndpixels * sizeof(double)); // mem for FLATS row to write
			   bpix = (double *) malloc(ndpixels * sizeof(double));  // mem for BIAS  row to write
			   cpix = (double *) malloc(ndpixels * sizeof(double));  // mem for DATA  row to write

			   if (apix == NULL || bpix == NULL || cpix == NULL) {
				   bail("Memory allocation error\n");
				}

			 } else {
					bail("Output file already exists %s\n",argv[4]);
			 }


			// Verify if the data image is the same size as the Flat and Bias. If not
			// check for a subrect keyword. Notice we check the 2D dimensions as it is assumed
			// that all images in a file are ofthe same size
			//
			if (( anaxes[0] != cnaxes[0] || anaxes[1] != cnaxes[1] )) {
				bail("Image Mismatch Size Error. Try using fitscopy to resize MasterBias or MasterFlat\n");
			}

			// loop over all rows of all of the data images
			for (firstpix[1] = 1; firstpix[1] <= (cnaxes[1]); firstpix[1]++) {

				bzero((void *) apix, ndpixels * sizeof(apix[0]));
				bzero((void *) bpix, ndpixels * sizeof(bpix[0]));
				bzero((void *) cpix, ndpixels * sizeof(cpix[0]));

				indexpix[1] = firstpix[1];

				if (fits_read_pix(mffptr, TDOUBLE, indexpix, ndpixels, NULL, bpix, NULL, &status)) {
							  bail("Failed to read Flat File row %ld \n",indexpix[1]);
				}

				if (fits_read_pix(mbfptr, TDOUBLE, indexpix, ndpixels, NULL, cpix, NULL, &status)) {
							  bail("Failed to read Bias File row %ld \n",indexpix[1]);
				}

				// This code will loop through each of the images in the data file and process the current row
				// this is done as it saves time reading in the master bias and flat for each row of each image
				for (firstpix[2] = 1; firstpix[2] <= cnaxes[2]; firstpix[2]++) {
				     if (fits_read_pix(datafptr, TDOUBLE, firstpix, ndpixels, NULL, apix, NULL, &status)) {
					     bail("Failed to read Flat File row %ld \n",firstpix[1]);
				     }

				    for(ii=0; ii< ndpixels; ii++) {
				            cval = (long) (apix[ii]-cpix[ii])/bpix[ii];
					    apix[ii] = (apix[ii]-cpix[ii])/bpix[ii];
				    }
				    fits_write_pix(outfptr, TDOUBLE, firstpix, ndpixels, apix, &status);
				}

			}


			// Read the DATE keyword value from the HDU
			ffgkey(outfptr, "DATE", keyval, comment ,&status);
			if (status) {
			   fits_report_error(stderr, status); /* print error message */
			   return(status);
			}
			length = strlen(keyval);

			  to=strndup(keyval+1, length-2);


			/* Converts a character buffer to correct numbers */
			ffs2tm(to, &year, &month, &day, &hour, &minute, &second, &status);

			/* Report error if status not = 0 */
			if (status) {
			   fits_report_error(stderr, status); /* print error message */
			   return(status);
			}

			/* convert date to julian date */
			jd = pr_julian_date(year,month,day,hour,minute,second);
			/* write the date to the output file */

			// Update the new file DATE keyword value with the Julian Date
			if (pr_update_date(outfptr, jd, &status) > 0)
			{
				return (status);
			}

            	fits_close_file(outfptr, &status);
			// Close the input data file
    		fits_close_file(datafptr, &status);
			if (status) {
           		fits_report_error(stderr, status); // print error message
           		bail(NULL);
    		}

    		free(apix);
		free(bpix);
            	free(cpix);

    }


    // Close all of the files

    fits_close_file(mffptr,  &status);
    fits_close_file(mbfptr,  &status);

	if (status) {
           fits_report_error(stderr, status); // print error message
           bail(NULL);
    }

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

/*-----------------------------------------------------------------
 * pr_update_date:
 * Parameters - fitsfile to update
 *              julian date to write to the DATE Keyword Value
 *              Status to indicate if operation successful
 *              If the Keyword already exists update the value
 */
int pr_update_date( fitsfile *fptr, double jd, int *status)
{
    int timeref;
    char date[30], card[FLEN_CARD];
    char DString[30];
	sprintf(DString,"%lf",jd); /* conver double to string */

    if (*status > 0)           /* inherit input status value if > 0 */
        return(*status);

    ffpky(fptr, TDOUBLE, "JD",&jd, "Original file creation date as Julian Date", status); /* Update the HDU */

    return(*status);
}

/*-----------------------------------------------------------------
 *	pr_julian_date:
 *  			computes the julian decimal date (j_date) from
 *				the gregorian calendar date.
 *				The Gregorian calendar reform occurred
 *				on 15 Oct. 1582.  This is 05 Oct 1582 by the
 *				julian calendar.
 *	Input:  Day, Month, Year, Hour, Minute, and second.
 *
 *	Output: the j_date will be the return value of the function.
 *
 *	Reference: Astronomial formulae for calculators, meeus, p 23
 *	from fortran program by F. Espenak - April 1982 Page 276,
 *	50 Year canon of solar eclipses: 1986-2035
 */

double pr_julian_date (int year,int month,int day,  int hour, int minute, double second)
{

	/* decimal day fraction	*/

	double	j_date; /* julian decimal date, 0 = 01 Jan -4712 12 HR UT */
	double frac, gyr;
	int jdn;

	int a,y,m;

    a = ((14-month)/12);
    y = year +4800 -a;
    m = month + (12*a)-3;
    jdn = day + (((153*m)+2)/5) + (365*y) + (y/4) - (y/100) + (y/400) -32045;
    frac = (((double)hour - 12)/24.0) + ((double)minute/1440) + (second/86400);

    j_date = (double) jdn + frac;
	return (j_date);
}
