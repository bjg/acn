#include <string.h>
#include "fitsio.h"
#include <sys/types.h>
#include <sys/dir.h>
#include <sys/param.h>
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <sys/resource.h>
#include <math.h>

#include <strings.h>

int     centroid(double *x, double *y, double *subrectarray, int xpos, int ypos,int boxdims,int threshold);
int     calc_magnitude(double centx, double centy, double *subrectarray, int xpos, int ypos,int boxdims, double radius, double skyB);
int 	compare_doubles (const void *X, const void *Y);

int     skybackground(double centx, double centy, double *subrectarray, int xpos, int ypos,int boxdims, double radius, double *median);
double  euclidian_dist(int pixelxpos,int pixelypos,double centx,double centy);
extern  int alphasort();

double C = 24;   // This is a key constant used in the Magnitude calculation

/*
** centroid:find a centroid from a rough x,y co-ordinate.
*		Paul Doyle 2012, Dublin Institute of Technology
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
    fprintf(stderr, "Usage: gen-mag directory x y boxsize \n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples: \n");
    fprintf(stderr, "  centroid ./objectfiledir 100 101 50 \n");
}

int main(int argc, char *argv[])
{
    fitsfile *datafptr;  /* FITS file pointers */
    int status = 0;      /* CFITSIO status value MUST be initialized to zero! */
	int anaxis, ii, Threshold = 660;
    
    long npixels = 1, ndpixels=1, fpixel[3] = {1,1,1},lpixel[3] = {1,1,1},inc[2]={0,1}, anaxes[3] = {1,1,1};    
    
    double *apix,radius=0; 
    int file_select();

    // variables to help read list of files
    int count,i,path_max = pathconf(".", _PC_NAME_MAX);
    struct direct **files;
    char fullfilename[path_max];  //to store path and filename
    double By=0,Bx=0,skyb=0;


    // Verify we have the correct number of parameters
    if (argc != 5)       { usage();  exit(0); }
    if (argv[1] == NULL) { usage();bail("Error - need x & y co-ordinates \n");  }
    if (argv[2] == NULL) { usage();bail("Error in arguments x");  }
    if (argv[3] == NULL) { usage();bail("Error in arguments y");  }
    if (argv[4] == NULL) { usage();bail("Error in box size ");  }

    int xpos = atoi(argv[2]); //convert the x input value to an integer
    int ypos = atoi(argv[3]); //convert the y input value to an integer
    int boxsize = atoi(argv[4]); //convert the boxsize input value to an integer

    count =  scandir(argv[1], &files, file_select, alphasort);
    printf ("Processing %d files \n",count);

	//
	// Process each of the data files, Remembering that the data files may be Cubed files
	//
 //       for (i=0;i<count;++i) {
    
           for (i=0;i<1;++i) {

			// Open the input files
        	snprintf(fullfilename, path_max - 1, "%s%s", argv[1], files[i]->d_name);
            printf ("Processing File...%s\n",fullfilename);
            fits_open_file(&datafptr, fullfilename, READONLY, &status); // open input images
        	if (status) {
        	   fits_report_error(stderr, status); // print error message
        	   bail(NULL);
        	}

			// Check the dimension of the DATA file */
			// cnaxis give the dimensions */
			fits_get_img_dim(datafptr, &anaxis, &status);  // read dimensions

			// Next we get the dimension filled in our 3D array anaxes
			fits_get_img_size(datafptr, 3, anaxes, &status);
			if (status) {
			   fits_report_error(stderr, status); /* print error message */
			   return(status);
			}

            ndpixels = boxsize*boxsize;  // 50 rows and 50 columns this is the number of pixels to store.
            apix = (double *) malloc(ndpixels * sizeof(double)); // mem for rectangle dataset

            if (apix == NULL) {
				   bail("Memory allocation error\n");
            }

				bzero((void *) apix, ndpixels * sizeof(apix[0]));

            fpixel[0] = xpos-boxsize/2;     // set up coordinates for a subrect which is 
            fpixel[1] = ypos-boxsize/2;     // 50 x50 width and height around the selected x/y coordinate provided
            lpixel[0] = xpos+boxsize/2-1;   // need to put in code to verify they box size is OK
            lpixel[1] = ypos+boxsize/2-1;
            inc[0] = inc[1] = 1;    // read all data pixels, don't skip any
                
            if (fpixel[0] < 1 || fpixel[1] < 1 || lpixel[0] < 1 || lpixel[1] < 1) 
                bail ("Not able to get a box area around the x,y coordinate");
                                
            // This code will loop through each of the images in the data Cube and process the current subrect identified
            for (fpixel[2] = 1; fpixel[2] <= anaxes[2]; fpixel[2]++) {
                printf ("Working on Image %ld \n",fpixel[2]);
                if (fits_read_subset(datafptr,TDOUBLE,fpixel,lpixel,inc,NULL,apix,NULL, &status)) {
                    fits_report_error(stderr, status); // print error message
                    bail("Failed to read subset of the image \n");
                }
                printf ("Radius     X        Y          S       I       SkyB       Flux \n");
                centroid (&Bx,&By,apix,xpos, ypos,boxsize,Threshold);
                //printf ("X,Y = %f %f \n",Bx,By);  // <-- Use this line to show centroid values
                skybackground(Bx,By,apix,xpos,ypos,boxsize,radius,&skyb);
                //Generate software aperature of varying sizes
                //
                for (radius = 1;radius<18;radius++) {
                    //printf ("radius = %f | ",radius);
                    calc_magnitude(Bx,By,apix,xpos,ypos,boxsize,radius,skyb);
                }
            }
			// Close the input data file
    		fits_close_file(datafptr, &status);
			if (status) {
           		fits_report_error(stderr, status); // print error message
           		bail(NULL);
    		}
    		free(apix);
    }

    exit(0);
}



/*
    skybackground: This function will create an annulus and a dannulus around the centre of the object
                   and calculate the skybackground by finding the MEDIAN value all all pixels found (exludes
                   partial pixels)
 */
int skybackground(double centx, double centy, double *subrectarray, int xpos, int ypos,int boxdims, double radius,double *median ) {

    float By = 0;
    float Bx = 0;
    float dist = 0;
    int pixelmaskcounter = 0;
    int ii,b,pixelxpos=0,pixelypos=0,mask=0;
    int rowcounter,colcounter,Npix=0;
    float readval=0,annulus = 0, dannulus=0;
    double *bpix,I = 0, S = 0,skyB = 50,Magnitude;
    double *dpix,medianval=0;

    annulus = radius+10;
    dannulus = annulus+15;
    
    if (dannulus > boxdims/2) {
        printf ("radius = %f, annulus = %f, dannunlus = %f, boxdims/2 = %d \n",radius, annulus,dannulus,boxdims/2);
        bail ("Need larger box around centre point to compute dannulus\n");
    }
    bpix = (double *) malloc(boxdims*boxdims * sizeof(double));  // mem for Mask
    dpix = (double *) malloc(boxdims*boxdims * sizeof(double));  // max amount of pixels requried
   
    if (bpix == NULL || dpix == NULL) 
        bail("Memory allocation error\n");

    
    
    //printf ("\n");
    for(ii=0; ii< boxdims; ii++) {   // loop over the rows
        rowcounter=0;
        for ( b=0;b<boxdims;b++) {   // loop over elements in the rows

            
            pixelxpos = b + xpos-boxdims/2;  // Get the x point in the frame not just the subrect
            pixelypos = ii + ypos-boxdims/2; // Get the y point in the frame not just the subrect
            
            dist = euclidian_dist(pixelxpos,pixelypos,centx,centy);  // Find dist to that pixel from the centre
            
            //
            // This section calculates the % of the pixel intensity to use.
            //
            if (dist < (dannulus -0.5) && (dist > (annulus +0.5)))  {
                //printf ("%f\n",subrectarray[ii*boxdims+b]);                   // inside the dannulus - full pizel
                            // Keep track of the numebr of pixels in the dannulus
                dpix[Npix++] = subrectarray[ii*boxdims+b];  // store the pixel value for sorting later
            }
            
            //bpix[ii*boxdims+b] = mask;      // We can use this mask to multiple by pixel intensity 
            
            // printf ("%d,%d -->%f<-- mask %d <<",pixelxpos,pixelypos,dist,mask);
        }
    }
        
    qsort(dpix,Npix,sizeof(double),compare_doubles);  // Sort all of the values
        
    // Find the median value
    if ( Npix % 2 == 0 )
            *median  = (dpix[(Npix/2)-1] + dpix[(Npix/2)])/2; // Even
    else
            *median  = (dpix[((Npix+1)/2) -1]); // Odd
   
        
        // printf ("Median = %f \n",*median);    
        
    free(bpix);
    free(dpix);
    
    return 0; // return value when all OK.
    
    
    
}


/*
 
 */

int calc_magnitude(double centx, double centy, double *subrectarray, int xpos, int ypos,int boxdims, double radius,double skyB) {

    float By = 0;
    float Bx = 0;
    float dist = 0;
    int pixelmaskcounter = 0;
    int ii,b,pixelxpos=0,pixelypos=0,mask=0;
    int rowcounter,colcounter;
    float readval=0;
    double *bpix,Npix=0,I = 0, S = 0,Magnitude;
    
    bpix = (double *) malloc(boxdims*boxdims * sizeof(double)); // mem for Mask
    //printf ("\n");
    for(ii=0; ii< boxdims; ii++) {   // loop over the rows
        rowcounter=0;
        for ( b=0;b<boxdims;b++) {   // loop over elements in the rows
            // Generate the bitmask using a threshold value of 660
            pixelxpos = b + xpos-boxdims/2;
            pixelypos = ii + ypos-boxdims/2;
            
            dist = euclidian_dist(pixelxpos,pixelypos,centx,centy);
            
            //
            // This section calculates the % of the pixel intensity to use.
            //
            if (dist == radius) 
                mask = 0.5;                 // On the aperature use 1/2 the pixel values
            else if (dist < (radius -0.5))  
                mask = 1;                   // inside the aperture use all the pixel values
            else if (dist > (radius + 0.5)) 
                mask = 0;                   // outside the aperture don't use pixel values
            else 
                mask = radius+0.5-dist;     // Calculate the partial pixel % further away is smaller
                
            //bpix[ii*boxdims+b] = mask;      // We can use this mask to multiple by pixel intensity 
            Npix += mask;                   // Keep track of the numebr of pixels in the aperture
            
            
            S+=subrectarray[ii*boxdims+b]*mask;  // Sum of the values within the aperature. 
            
           // printf ("%d,%d -->%f<-- mask %d <<",pixelxpos,pixelypos,dist,mask);
        }
    }
    I = S - (skyB * Npix);
    Magnitude = (-2.5 * log10(I)) + C;   
    printf ("%7.0f %7.4f %7.4f %7.4f %7.4f %7.2f %7.5f \n",radius,centx,centy,S,I,skyB, Magnitude);

    free(bpix);
    
    return 0; // return value when all OK.
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

int centroid(double *x, double *y, double *subrectarray, int xpos, int ypos,int boxdims, int threshold) {
    
    float By = 0;
    float Bx = 0;
    int pixelmaskcounter = 0;
    int ii,b;
    int rowcounter,colcounter;
    float readval=0;
    double *bpix;
  
    bpix = (double *) malloc(boxdims*boxdims * sizeof(double)); // mem for Mask

    for(ii=0; ii< boxdims; ii++) {
        rowcounter=0;
        for ( b=0;b<boxdims;b++) {
            // Generate the bitmask using a threshold value of 660
            readval = subrectarray[ii*boxdims+b];
            if (readval > threshold) {
                bpix[ii*boxdims+b] = 1;
                pixelmaskcounter++;
                rowcounter++;  // counter the number of events in the row
            }
            else
                bpix[ii*boxdims+b] = 0;
        }
        By += rowcounter * ii;   // By count of pixel row elements
    }
    
    // Process the colums to get the X centre point
    for(ii=0; ii< boxdims; ii++) {
        colcounter=0;
        for ( b=0;b<boxdims;b++) {
            // Read the bitmask 
            if (bpix[ii+(boxdims*b)] == 1)
                colcounter++;  // counter the number of events in the col
        }
        Bx += colcounter * ii;   // By count of pixel row elements
    }
    
    *x = (Bx/pixelmaskcounter)+xpos-(boxdims/2);   // return the global location in the frame 
    *y = (By/pixelmaskcounter)+ypos-(boxdims/2);
  
    free(bpix);
    return *x,*y;
    
}

// Caclulate the distance between a pixel point and the centre of the point source/star
double euclidian_dist(int pixelxpos,int pixelypos,double centx,double centy) {
    return sqrt((((double)pixelxpos - centx)*((double)pixelxpos - centx)) + (((double)pixelypos - centy)*((double)pixelypos - centy)));
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
    


int file_select(struct direct   *entry) {
                             
    char *ptr;
                             
    if ((strcmp(entry->d_name, ".")== 0) ||    (strcmp(entry->d_name, "..") == 0))
        return (FALSE);
                             
    /* Check for filename extensions */
    ptr = strrchr(entry->d_name, '.');
    return ((ptr != NULL) && (strcmp(ptr, ".fits") == 0));
}
