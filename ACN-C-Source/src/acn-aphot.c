#include <string.h>
#include <fitsio.h>
#include <sys/types.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/resource.h>
#include <math.h>
#include <errno.h>
#include <strings.h>

#define BUFSIZE 2056
#define MINRADIUS 2
#define MAXSTARS 15

#define C 24

/**
*
*  Function Prototypes
*
* 
* 
**/
static FILE *open_result_file(const char *prefix);
int centroid(double *x, double *y, double *subrectarray, double *mfarray,
	     double *mbarray, int xpos, int ypos, int boxdims,
	     int threshold);
int calc_magnitude(double centx, double centy, double *subrectarray,
		   double *mfarray, double *mbarray, int xpos, int ypos,
		   int boxdims, double radius, double skyB);
int compare_doubles(const void *X, const void *Y);
int skybackground(double centx, double centy, double *subrectarray,
		  double *mfarray, double *mbarray, int xpos, int ypos,
		  int boxdims, double annulusval, double dannulusval,
		  double radius, double *median);
double euclidian_dist(int pixelxpos, int pixelypos, double centx,
		      double centy);
extern int alphasort();
double xguessarray[MAXSTARS], yguessarray[MAXSTARS], radiusarray[MAXSTARS],
    annulusarray[MAXSTARS], dannulusarray[MAXSTARS], boxarray[MAXSTARS],
    thresholdarray[MAXSTARS];
int debug = 0;
char buf[BUFSIZE], *p, *token[7];
int cleanmode = 0, i = 0, sc = 0;
FILE *fp = NULL;		// This is used within multiple functions

/*
*      acn-aphot: Estimate of point sources magnitude values using either a cleaned data file, or raw 
*                 file which can be cleaned during the process. 
*
*        Paul Doyle 2012, Dublin Institute of Technology
*/

// 
// Exit the program providing an error message to the stderr
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


//
// Provide users with details on how to use the programme
//
void usage(void)
{
    fprintf(stderr,
	    "Usage: acn-aphot ./directory [-c ./masterflat ./masterbias] < ./config \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  acn-aphot ./objectfiledir < ./config\n");
    fprintf(stderr,
	    "  acn-aphot ./objectfiledir -c ./masterflat ./masterbias < ./config\n\n");

}

int main(int argc, char *argv[])
{
    fitsfile *datafptr, *mffptr, *mbfptr;	/* FITS file pointers */
    //char    buf[BUFSIZE]; // *p;   

    int status = 0;		/* CFITSIO status value MUST be initialized to zero! */
    int anaxis, bnaxis, cnaxis, j, ii, Threshold = 660;
    long npixels = 1, ndpixels = 1, fpixel[3] = { 1, 1, 1 }, lpixel[3] = {
    1, 1, 1}, inc[2] = {
    0, 1}, anaxes[3] = {
    1, 1, 1}, bnaxes[3] = {
    1, 1, 1}, cnaxes[3] = {
    1, 1, 1};

    double *apix, radius = 0;
    double *bpix[MAXSTARS];	// An Array of pointers for bias 
    double *cpix[MAXSTARS];	// An Array of pointers for flats

    int file_select();

    // variables to help read list of files
    int count, i = 0, path_max = pathconf(".", _PC_NAME_MAX);
    struct direct **files;
    char fullfilename[path_max];	//to store path and filename
    double By = 0, Bx = 0, skyb = 0;


    //
    // Verify we have the correct number of parameters
    //

    if (argc == 5) {		// Verify that we selected  cleanmode.
	if (strcmp(argv[2], "-c") == 0)
	    cleanmode = 1;
	else {
	    usage();
	    bail("Invalid parameters\n");
	}
    } else if (argc != 2) {
	usage();
	bail("Invalid parameters\n");
    }
    //
    // Parse the configuration file
    //
    while (fgets(buf, BUFSIZE, stdin) != NULL) {
	if (buf[0] == '!') {
	    ;			// We just ignore the comment lines
	} else {
	    errno = 0;
	    xguessarray[sc] = atof(strtok(buf, " "));
	    if (xguessarray[sc] == 0)
		bail("invalid config file X value");
	    yguessarray[sc] = atof(strtok(NULL, " "));
	    if (yguessarray[sc] == 0)
		bail("invalid config file Y value");
	    radiusarray[sc] = atof(strtok(NULL, " "));
	    if (radiusarray[sc] == 0)
		bail("invalid config file Radius value");
	    annulusarray[sc] = atof(strtok(NULL, " "));
	    if (annulusarray[sc] == 0)
		bail("invalid config file annulus value");
	    dannulusarray[sc] = atof(strtok(NULL, " "));
	    if (dannulusarray[sc] == 0)
		bail("invalid config file dannulus value");
	    boxarray[sc] = atof(strtok(NULL, " "));
	    if (boxarray[sc] == 0)
		bail("invalid config file boxarry value");
	    thresholdarray[sc] = atof(strtok(NULL, " "));
	    if (thresholdarray[sc] == 0)
		bail("invalid config file threshold value");
	    if (errno != 0)
		bail("Invalid config file detected");
	    sc++;		// Increase the star counter
	    if (sc > MAXSTARS)
		bail("too many stars to proces in the config file\n");
	}
    }

    if (sc < 1) {
	usage();
	bail("Config file did not contain star data points\n");
    } else {
	printf("Found %d stars \n", sc);
    }

    // If we are going to clean the image then we have to first check the dimensions of the supplied
    // Master Bias and Master Flat. 
    if (cleanmode == 1) {	// This means we have to read the master flat and bias 

	fits_open_file(&mffptr, argv[3], READONLY, &status);	// open master flat file
	if (status) {
	    fits_report_error(stderr, status);	// print error message
	    bail(NULL);
	}

	fits_open_file(&mbfptr, argv[4], READONLY, &status);	// open master bias file
	if (status) {
	    fits_report_error(stderr, status);	// print error message
	    bail(NULL);
	}
	// Read the flat and bias files to establish the dimensions are identical.
	fits_get_img_dim(mffptr, &bnaxis, &status);	// read dimensions of the file
	fits_get_img_size(mffptr, 3, bnaxes, &status);
	if (status) {
	    fits_report_error(stderr, status);	// print error message
	    bail(NULL);
	}
	if (bnaxis > 2)
	    bail("Error: Master Flat File %s in an images with > 2 dimensions and is not supported\n", argv[3]);

	// Verify that the Master Bias and Master Flat are the same dimensions
	fits_get_img_dim(mbfptr, &cnaxis, &status);	// read dimensions of each file
	fits_get_img_size(mbfptr, 3, cnaxes, &status);

	if (status) {
	    fits_report_error(stderr, status);	// print error message
	    bail(NULL);
	}
	if (cnaxis > 2)
	    bail("Error: Master Bias File %s in an images with > 2 dimensions and is not supported\n", argv[4]);

	// Bias and Master files should be the same size.
	if ((bnaxes[0] != cnaxes[0] || bnaxes[1] != cnaxes[1]))
	    bail("Error: input images don't have same size\n");

	// Read the master bias and the master flat files extracting the required values for each star
	for (i = 0; i < sc; ++i) {
	    fpixel[0] = xguessarray[i] - boxarray[i] / 2;	// set up coordinates for a subrect which is 
	    fpixel[1] = yguessarray[i] - boxarray[i] / 2;	// 50 x50 width and height around the selected x/y coordinate provided
	    lpixel[0] = xguessarray[i] + boxarray[i] / 2 - 1;	// need to put in code to verify they box size is OK
	    lpixel[1] = yguessarray[i] + boxarray[i] / 2 - 1;
	    inc[0] = inc[1] = 1;	// read all data pixels, don't skip any

	    ndpixels = boxarray[i] * boxarray[i];	// 50 rows and 50 columns this is the number of pixels to store.
	    bpix[i] = (double *) malloc(ndpixels * sizeof(double));	// mem for rectangle dataset
	    cpix[i] = (double *) malloc(ndpixels * sizeof(double));	// mem for rectangle dataset

	    if (bpix[i] == NULL || cpix[i] == NULL) {
		bail("Memory allocation error\n");
	    }
	    bzero((void *) bpix[i], ndpixels * sizeof(double));
	    bzero((void *) cpix[i], ndpixels * sizeof(double));

	    // Get the MF subrect value
	    if (fits_read_subset
		(mffptr, TDOUBLE, fpixel, lpixel, inc, NULL, bpix[i], NULL,
		 &status)) {
		fits_report_error(stderr, status);	// print error message
		bail("Failed to read subset Master Flat of the image \n");
	    }
	    // Get the MB subrect value
	    if (fits_read_subset
		(mbfptr, TDOUBLE, fpixel, lpixel, inc, NULL, cpix[i], NULL,
		 &status)) {
		fits_report_error(stderr, status);	// print error message
		bail("Failed to read subset Master Bias of the image \n");
	    }

	}
    }

    count = scandir(argv[1], &files, file_select, alphasort);
    printf("Processing %d files \n", count);
    //
    // Process each of the data files, Remembering that the data files may be Cubed files
    //
    for (i = 0; i < count; ++i) {

	// Open the input files
	snprintf(fullfilename, path_max - 1, "%s%s", argv[1],
		 files[i]->d_name);

	printf("Processing File...%s\n", files[i]->d_name);
	fits_open_file(&datafptr, fullfilename, READONLY, &status);	// open input images
	if (status) {
	    fits_report_error(stderr, status);	// print error message
	    bail(NULL);
	}
	// Check the dimension of the DATA file */
	// cnaxis give the dimensions */
	fits_get_img_dim(datafptr, &anaxis, &status);	// read dimensions

	// Next we get the dimension filled in our 3D array anaxes
	fits_get_img_size(datafptr, 3, anaxes, &status);
	if (status) {
	    fits_report_error(stderr, status);	/* print error message */
	    return (status);
	}

	if (cleanmode == 1) {
	    // Check if the bias and the object file are the same dimensions
	    if ((anaxes[0] != bnaxes[0] || anaxes[1] != bnaxes[1]))
		bail("Error: input images don't have same size\n");
	}

	fp = open_result_file(files[i]->d_name);

	// Loop through each of the stars found in the configuraiton file
	// They may contain different box sizes, x,y guesses, radius ranges and thresholds
	//
	// The looping means that for a particular star we see the data results for that star across each of the
	// images in the cube, we then go on the next start which we read through the cube. 
	// Need to check how we want to group the output. 
	for (j = 0; j < sc; j++) {

	    fprintf(fp,
		    "Processing star number %d in configuration file\n\n",
		    j + 1);
	    ndpixels = boxarray[j] * boxarray[j];	// 50 rows and 50 columns this is the number of pixels to store.
	    apix = (double *) malloc(ndpixels * sizeof(double));	// mem for rectangle dataset

	    if (apix == NULL) {
		bail("Memory allocation error\n");
	    }
	    bzero((void *) apix, ndpixels * sizeof(apix[0]));

	    fpixel[0] = xguessarray[j] - boxarray[j] / 2;	// set up coordinates for a subrect which is 
	    fpixel[1] = yguessarray[j] - boxarray[j] / 2;	// 50 x50 width and height around the selected x/y coordinate provided
	    lpixel[0] = xguessarray[j] + boxarray[j] / 2 - 1;	// need to put in code to verify they box size is OK
	    lpixel[1] = yguessarray[j] + boxarray[j] / 2 - 1;
	    inc[0] = inc[1] = 1;	// read all data pixels, don't skip any

	    if (fpixel[0] < 1 || fpixel[1] < 1 || lpixel[0] < 1
		|| lpixel[1] < 1)
		bail("Not able to get a box area around the x,y coordinate %d %d %d %d\n", fpixel[0], fpixel[1], lpixel[0], lpixel[1]);

	    // This code will loop through each of the images in the data Cube and process the current subrect identified
	    for (fpixel[2] = 1; fpixel[2] <= anaxes[2]; fpixel[2]++) {

		fprintf(fp, "Working on Image %ld \n", fpixel[2]);

		if (fits_read_subset
		    (datafptr, TDOUBLE, fpixel, lpixel, inc, NULL, apix,
		     NULL, &status)) {
		    fits_report_error(stderr, status);	// print error message
		    bail("Failed to read subset of the image \n");
		}
		fprintf(fp,
			"Radius     X        Y          S       I       SkyB       Mag Estimate \n");

		centroid(&Bx, &By, apix, bpix[j], cpix[j], xguessarray[j],
			 yguessarray[j], boxarray[j], thresholdarray[j]);

		//Generate software aperature of varying sizes
		for (radius = MINRADIUS; radius < radiusarray[0]; radius++) {
		    skybackground(Bx, By, apix, bpix[j], cpix[j],
				  xguessarray[j], yguessarray[j],
				  boxarray[j], annulusarray[j],
				  dannulusarray[j], radius, &skyb);
		    // Sky background changes as radius moves and pushes out the annulus
		    calc_magnitude(Bx, By, apix, bpix[j], cpix[j],
				   xguessarray[j], yguessarray[j],
				   boxarray[j], radius, skyb);
		}
	    }
	    free(apix);
	}

	fclose(fp);

	// Close the input data file
	fits_close_file(datafptr, &status);
	if (status) {
	    fits_report_error(stderr, status);	// print error message
	    bail(NULL);
	}
    }
    if (cleanmode == 1) {
	for (i = 0; i < sc; i++) {
	    free(bpix[i]);
	    free(cpix[i]);
	}

	fits_close_file(mffptr, &status);
	fits_close_file(mbfptr, &status);
	if (status) {
	    fits_report_error(stderr, status);	// print error message
	    bail(NULL);
	}
    }
    for (i = 0; i < count; i++) {
	free(files[i]);
    }
    free(files);
    exit(0);
}

static FILE *open_result_file(const char *prefix)
{
    const char *suffix = ".result";
    char *filename = strdup(prefix);
    filename = realloc(filename, strlen(prefix) + strlen(suffix) + 1);
    strcat(filename, suffix);
    FILE *fp = fopen(filename, "w");
    free(filename);
    return fp;
}


/*
    skybackground: This function will create an annulus and a dannulus around the centre of the object
                   and calculate the skybackground by finding the MEDIAN value all all pixels found (exludes
                   partial pixels)
 */
int
skybackground(double centx, double centy, double *subrectarray,
	      double *mfarray, double *mbarray, int xpos, int ypos,
	      int boxdims, double annulusval, double dannulusval,
	      double radius, double *median)
{

    float By = 0;
    float Bx = 0;
    float dist = 0;
    int pixelmaskcounter = 0;
    int ii, b, pixelxpos = 0, pixelypos = 0, mask = 0;
    int rowcounter, colcounter, Npix = 0;
    float readval = 0, annulus = 0, dannulus = 0;
    double *epix, I = 0, S = 0, skyB = 50, Magnitude;
    double *dpix, medianval = 0;

    annulus = radius + annulusval;
    dannulus = annulus + dannulusval;

    if (dannulus > boxdims / 2) {
	printf
	    ("radius = %f, annulus = %f, dannunlus = %f, boxdims/2 = %d \n",
	     radius, annulus, dannulus, boxdims / 2);
	bail("Need larger box around centre point to compute dannulus\n");
    }
    epix = (double *) malloc(boxdims * boxdims * sizeof(double));	// mem for Mask
    dpix = (double *) malloc(boxdims * boxdims * sizeof(double));	// max amount of pixels requried

    if (epix == NULL || dpix == NULL)
	bail("Memory allocation error\n");

    for (ii = 0; ii < boxdims; ii++) {	// loop over the rows
	rowcounter = 0;
	for (b = 0; b < boxdims; b++) {	// loop over elements in the rows


	    pixelxpos = b + xpos - boxdims / 2;	// Get the x point in the frame not just the subrect
	    pixelypos = ii + ypos - boxdims / 2;	// Get the y point in the frame not just the subrect

	    dist = euclidian_dist(pixelxpos, pixelypos, centx, centy);	// Find dist to that pixel from the centre

	    //
	    // This section calculates the % of the pixel intensity to use.
	    //
	    if (dist < (dannulus - 0.5) && (dist > (annulus + 0.5))) {
		if (cleanmode == 1)
		    dpix[Npix++] = (subrectarray[ii * boxdims + b] - mbarray[ii * boxdims + b]) / mfarray[ii * boxdims + b];	// store the pixel value for sorting later
		else
		    dpix[Npix++] = subrectarray[ii * boxdims + b];	// store the pixel value for sorting later
	    }
	}
    }

    qsort(dpix, Npix, sizeof(double), compare_doubles);	// Sort all of the values

    // Find the median value
    if (Npix % 2 == 0)
	*median = (dpix[(Npix / 2) - 1] + dpix[(Npix / 2)]) / 2;	// Even
    else
	*median = (dpix[((Npix + 1) / 2) - 1]);	// Odd

    free(epix);
    free(dpix);

    return 0;			// return value when all OK.
}


/*
 
 */

int
calc_magnitude(double centx, double centy, double *subrectarray,
	       double *mfarray, double *mbarray, int xpos, int ypos,
	       int boxdims, double radius, double skyB)
{

    float By = 0;
    float Bx = 0;
    float dist = 0;
    int pixelmaskcounter = 0;
    int ii, b, pixelxpos = 0, pixelypos = 0, mask = 0;
    int rowcounter, colcounter;
    float readval = 0;
    double *bpix, Npix = 0, I = 0, S = 0, Magnitude;

    bpix = (double *) malloc(boxdims * boxdims * sizeof(double));	// mem for Mask
    for (ii = 0; ii < boxdims; ii++) {	// loop over the rows
	rowcounter = 0;
	for (b = 0; b < boxdims; b++) {	// loop over elements in the rows
	    // Generate the bitmask using a threshold value of 660
	    pixelxpos = b + xpos - boxdims / 2;
	    pixelypos = ii + ypos - boxdims / 2;
	    dist = euclidian_dist(pixelxpos, pixelypos, centx, centy);
	    //
	    // This section calculates the % of the pixel intensity to use.
	    //
	    if (dist == radius)
		mask = 0.5;	// On the aperature use 1/2 the pixel values
	    else if (dist < (radius - 0.5))
		mask = 1;	// inside the aperture use all the pixel values
	    else if (dist > (radius + 0.5))
		mask = 0;	// outside the aperture don't use pixel values
	    else
		mask = radius + 0.5 - dist;	// Calculate the partial pixel % further away is smaller

	    Npix += mask;	// Keep track of the numebr of pixels in the aperture

	    if (cleanmode == 1)
		S += ((subrectarray[ii * boxdims + b] -
		       mbarray[ii * boxdims + b]) / mfarray[ii * boxdims +
							    b]) * mask;
	    else
		S += subrectarray[ii * boxdims + b] * mask;
	}
    }
    I = S - (skyB * Npix);
    Magnitude = (-2.5 * log10(I)) + C;
    fprintf(fp, "%7.0f %7.4f %7.4f %7.6f %7.6f %7.2f %7.5f \n", radius,
	    centx, centy, S, I, skyB, Magnitude);

    free(bpix);

    return 0;			// return value when all OK.
}


//
// Given an array of values representing a rectangular part of the image containing 
// a point source find the X & Y centroid value
// Need to pass in the following parameters 
// 
//  x,y - these are used to pass the centroid values back to the calling program
//  subrect - This is the array returned from the fits_read_subset function containing values from a specific region
//  xpos,ypos - this is the offset X,Y positions use to calculate where the real centroid value is. This was our initial guess
//  boxdims - width & heigh of the box (assumed to be a square - this is very important!
//  threshold - value used to determine what value is cutoff for use in the mask.
//  
//
//  Paul Doyle -  March 2012

int
centroid(double *x, double *y, double *subrectarray, double *mfarray,
	 double *mbarray, int xpos, int ypos, int boxdims, int threshold)
{

    float By = 0;
    float Bx = 0;
    int pixelmaskcounter = 0;
    int ii, b;
    int rowcounter, colcounter;
    float readval = 0;
    double *bpix;

    bpix = (double *) malloc(boxdims * boxdims * sizeof(double));	// mem for Mask

    for (ii = 0; ii < boxdims; ii++) {
	rowcounter = 0;
	for (b = 0; b < boxdims; b++) {
	    // Generate the bitmask using a threshold value of 660
	    if (cleanmode == 1)
		readval =
		    (subrectarray[ii * boxdims + b] -
		     mbarray[ii * boxdims + b]) / mfarray[ii * boxdims +
							  b];
	    else
		readval = subrectarray[ii * boxdims + b];

	    if (readval > threshold) {
		bpix[ii * boxdims + b] = 1;
		pixelmaskcounter++;
		rowcounter++;	// counter the number of events in the row
	    } else
		bpix[ii * boxdims + b] = 0;
	}
	By += rowcounter * ii;	// By count of pixel row elements
    }

    // Process the colums to get the X centre point
    for (ii = 0; ii < boxdims; ii++) {
	colcounter = 0;
	for (b = 0; b < boxdims; b++) {
	    // Read the bitmask 
	    if (bpix[ii + (boxdims * b)] == 1)
		colcounter++;	// counter the number of events in the col
	}
	Bx += colcounter * ii;	// By count of pixel row elements
    }

    *x = (Bx / pixelmaskcounter) + xpos - (boxdims / 2);	// return the global location in the frame 
    *y = (By / pixelmaskcounter) + ypos - (boxdims / 2);

    free(bpix);

    return *x, *y;

}

// Caclulate the distance between a pixel point and the centre of the point source/star
double
euclidian_dist(int pixelxpos, int pixelypos, double centx, double centy)
{
    return
	sqrt((((double) pixelxpos - centx) * ((double) pixelxpos -
					      centx)) +
	     (((double) pixelypos - centy) * ((double) pixelypos -
					      centy)));
}

int compare_doubles(const void *X, const void *Y)
{
    double x = *((double *) X);
    double y = *((double *) Y);

    if (x > y) {
	return 1;
    } else {
	if (x < y) {
	    return -1;
	} else {
	    return 0;
	}
    }
}


int file_select(struct direct *entry)
{

    char *ptr;

    if ((strcmp(entry->d_name, ".") == 0)
	|| (strcmp(entry->d_name, "..") == 0))
	return (FALSE);

    /* Check for filename extensions */
    ptr = strrchr(entry->d_name, '.');
    return ((ptr != NULL) && (strcmp(ptr, ".fits") == 0));
}
