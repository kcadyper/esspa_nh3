#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include "mc_custom.h"

#define RAND_MAX 2147483647

#define FALSE 0
#define TRUE 1

#ifndef NULL
#define NULL 0
#endif

#define MINTAU 0.01
#define PI 3.14159265
#define DEG2RAD (PI/180.0)
#define TCOSMIC 2.7
#define BIG 1.0e6


#define IBMAX (IBOXES-1)
#define JBMAX (JBOXES-1)
#define KBMAX (KBOXES-1)

#define IWMAX IBOXES
#define JWMAX JBOXES
#define KWMAX KBOXES

#define IPIXELS (IBOXES*PPB)
#define JPIXELS (JBOXES*PPB)

#define DPX (DSH/PPB)      /* dimensions of image pixels */
#define DPY (DSH/PPB)


/* viewing angle */
extern float phi0, theta0 ;

/* initial propagation unit vector */
extern float k10, k20, k30 ;

/* illumination unit vector for visible simulations*/
extern float k1_illum, k2_illum, k3_illum ;

/* current propagation unit vector */
extern float k1, k2, k3 ;

/* current position vector */
extern float d1, d2, d3 ;

/* current box */
extern int ib,jb,kb ;

/* cloud-free TB for calculation of P, S */
extern float tv0, th0 ;

/* upper and lower brightness temperature limits for display purposes */
extern float tmin,tmax ;
extern float pmin,pmax,smin,smax  ;

/* number of photons to use per image pixel */
extern int nmax ;

/* reflectivity experienced by photon during 
most recent interaction with surface */
extern float reflecth, reflectv ;

/* single-scattering albedo experienced by photon during 
most recent scattering event */
extern float ssa ;

/* brightness temperature images for vertical and horizontal pol. */

extern float tbh[IPIXELS][JPIXELS] , tbv[IPIXELS][JPIXELS] ;

/* flags indicating locations of valid data in image */
extern int pixflag[IPIXELS][JPIXELS] ;

/* photon count for each pixel */
extern int ncount[IPIXELS][JPIXELS] ;

/* volume extinction coefficient */
extern float kext[IBOXES][JBOXES][KBOXES] ;

/* gridbox coordinates */
extern float x[IBOXES+1] , y[JBOXES+1] , z[KBOXES+1] ;

/* box into which photon will exit if it continues on current path */
extern int ibnew,jbnew,kbnew ;

/* coordinates of photon at entry to next box */
extern float d1new,d2new,d3new ;

/* cumulative temperature of photon */
extern float tbphoth, tbphotv ;

/* cumulative transmittance of photon */
extern float tauv, tauh ;

/* optical distance to next scattering */
extern float sigma ; 

/* current image pixel position */
extern int ipix, jpix ;

/* geometric distance to exit */
extern float s ;

/* whether photon has been scattered */
extern int scattered ;

/* domain averaged brightness temperature results */
extern float tbvave,tbhave,pave,save;

/* microwave frequency used in simulation */
/*extern float freq ;*/
static float freq ;

extern double local_temp() ;
extern double surface_temp() ;
extern double uran() ;
extern double get_rand_pathlength() ;
extern int par_init() ;
extern int main_loop(int);

/* input filename */
extern char *infile ;


