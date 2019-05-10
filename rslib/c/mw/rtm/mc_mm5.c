#define _POSIX_SOURCE
#include <stdio.h>
#include <fcntl.h>
#include <math.h>
#include "mc.h"

/*static float freq ;*/

static char outfile[256] ;


float Temp[IBOXES][JBOXES][KBOXES] ;
float ssag[IBOXES][JBOXES][KBOXES] ;
float assym[IBOXES][JBOXES][KBOXES] ;

int include_cloud ;
int include_rain ;
int include_graupel;
int include_agg ;
int include_ice ;
int include_gas ;
int include_vapor ;

int aa;

double local_temp()
{
  return((double) Temp[ib][jb][kb]) ;
}

double surface_temp()
{
  double t ;
  /* use lower atmospheric temp as surface temp too */
  t = Temp[ib][jb][0] ;
  return(t) ;
}


get_scattering_pars()
{
  float phi, costhet, sinthet ;
  float sinthetprim ;
  float sinphiprim,cosphiprim ;
  float phiprim,den ;
  float v1,v2,v3 ;
  float h1,h2,h3 ;
  float x ;

  float g ;

  static float gold = -10.0 ;
  static float t1,t2,t4,t5 ;
  float kk1,kk2,kk3 ;

  double r,t3,t6,costhetprim ;

  ssa = ssag[ib][jb][kb] ;


#ifdef ISOTROPIC
  phi = uran()*2.0*PI ;
  costhet = 2.0*uran() - 1.0 ;
  if (costhet >  1.0) costhet =  1.0 ;
  if (costhet < -1.0) costhet = -1.0;
  sinthet = sqrt(1.0 - costhet*costhet) ;

  k1 = cos(phi)*sinthet ;
  k2 = sin(phi)*sinthet ;
  k3 = costhet ;
#endif

#ifdef HENYEYGREENSTEIN

  /* random number between 0.0 and 1.0 */
  r = uran() ;

  g = assym[ib][jb][kb] ;

  t4 = 1.0 - g ;
  t1 = t4*t4 ;
  t2 = g + g ;
  t5 = t1 + t2  ;

  t6 = g*r ;
  t3 = t4 + t6 ;
  costhetprim = (2.0*r*t5*t3 - t1)/(4.0*t6*t3 + t1) ;
  if (costhetprim > 1.0) costhetprim = 1.0 ;
  if (costhetprim< -1.0) costhetprim = -1.0 ;

  phiprim = uran()*2.0*PI ;

  x = 1.0-costhetprim*costhetprim ;
  if (x < 0.0) x = 0.0 ;
  sinthetprim = sqrt(x) ;
  cosphiprim = cos(phiprim) ;
  sinphiprim = sin(phiprim) ;


  den = sqrt(k1*k1 + k2*k2) ;
  if (den <= 1.0e-6 ) {
    h1 = 0.0 ;
    h2 = 1.0 ;
    h3 = 0.0 ;
    v1 = 1.0 ;
    v2 = 0.0 ;
    v3 = 0.0 ;
  } else {
    h1 = k2/den ;
    h2 = -k1/den ;
    h3 = 0.0 ; 
    v1 = h2*k3 ;
    v2 = -h1*k3 ;
    v3 = h1*k2-h2*k1 ;
    den = sqrt(v1*v1 + v2*v2 + v3*v3) ;
    v1 /= den ;
    v2 /= den ;
    v3 /= den ;
  }

  k1 = k1*costhetprim + sinthetprim*(v1*cosphiprim + h1*sinphiprim) ;
  k2 = k2*costhetprim + sinthetprim*(v2*cosphiprim + h2*sinphiprim) ;
  k3 = k3*costhetprim + sinthetprim*(v3*cosphiprim + h3*sinphiprim) ;
  den = sqrt(k1*k1 + k2*k2 + k3*k3) ;
  k1 /= den ;
  k2 /= den ;
  k3 /= den ;
#endif

}

get_reflection_pars()
{

  float mu,phi,sinthet ;
  float f,tk,theta,ssw,ev,eh;

#ifdef SPECULAR

  k3 = -k3 ;  

  /* assume Fresnel emissivity from ocean surface */

  f = freq ;
  tk = surface_temp() ;
  ssw = 36.0 ;
  theta = acos(k3)/DEG2RAD ;
  spemiss_(&freq,&tk,&theta,&ssw,&ev,&eh) ;


#endif
#ifdef LAMBERTIAN

  phi = uran()*2.0*PI ;   
  mu = 2.0*uran() - 1.0 ;
  if (mu >  1.0) mu =  1.0 ;
  if (mu < -1.0) mu = -1.0;
  sinthet = sqrt(1.0 - mu*mu) ;
  k1 = cos(phi)*sinthet ;
  k2 = sin(phi)*sinthet ;
  k3 = mu ;
  if (k3 < 0) k3 = -k3 ;   

  /* ad hoc emissivity values for land surface */
  ev = 0.95 ;
  eh = 0.85 ;

#endif

  if (scattered) {       /* if photon has been scattered at least once */
    ev = 0.5*(ev + eh) ;  /* then polarization is no longer known */
    eh = ev ;
  }
  reflectv = 1.0 - ev ;
  reflecth = 1.0 - eh ;
}



init_optics()
{

  float absair, abswv, abscld ;
  float delz ;
  float kext_total,ksca_total,g_total,ss ;
  int i,j,k, ix,jx,kx, it,jt,kt ;
  int ifreq ;
  float  rext,rssa,rg, 
    gext,gssa,gg, 
    sext,sssa,sg, 
    iext,issa,ig ;

  FILE *tfile ;
  char cdum[256], timestring[256], string[256] ;
  float dsh ;

  /* infile name was previously passed as argument to main() */

  tfile = fopen(infile, "r") ;
  if (tfile == (FILE *) NULL) {
    printf("Couldn't open file %s \n",infile) ;
    exit(-1) ;
  }

  sscanf(infile,"%*9c%6d",&ifreq ) ;
  freq = ifreq/1000.0 ;
  /*  printf("Freq. = %.3f GHz \n",freq) ; */
  
  fgets(string,255,tfile) ;
  /*printf("%s %45c \n",infile,string);*/

  /*  printf("Reading file %s containing %45c \n",infile,string) ;*/
  fgets(string,255,tfile) ;
  sscanf(string,"%8s",timestring) ;
  /*  printf("Date/time = %8c \n",timestring) ; */

  /* prompt for hydrometeor types to include */

  printf("Include Gaseous Abs.? (0,1): ") ;
  scanf("%d",&include_gas) ;
  if (include_gas) {
    printf("Include Vapor? (0,1): ") ;
    scanf("%d",&include_vapor) ;
  } else {
    include_vapor = FALSE ;
  }
  printf("Include Cloud? (0,1): ") ;
  scanf("%d",&include_cloud) ;
  printf("Include Rain? (0,1): ") ;
  scanf("%d",&include_rain) ;
  printf("Include Graupel? (0,1): ") ;
  scanf("%d",&include_graupel) ;
  printf("Include Snow? (0,1): ") ;
  scanf("%d",&include_agg) ;
  printf("Include Ice? (0,1): ") ;
  scanf("%d",&include_ice) ;


  /* create informative output file name */

  sprintf(outfile,"%s_%1d%1d%1d%1d%1d%1d%1d.tbs",infile,
	  include_gas,include_vapor,include_cloud,include_rain,include_graupel,include_agg,include_ice) ;
  printf("Output file name = %s\n", outfile) ;


  /* begin reading input grids */

  fgets(string,255,tfile) ;
  sscanf(string,"%*15c %d %d %d ",&ix,&jx,&kx) ;
  printf("Grid dimensions: %d %d %d \n",ix,jx,kx) ;
  
  if (ix != IBOXES || jx != JBOXES || kx != KBOXES) {
    printf("Input grid dimensions not compatible with values in mc_custom.h.\n Edit latter file and recompile.\n") ;
    exit(-1) ;
  }

  fgets(string,255,tfile) ;
  sscanf(string,"%*21c %f",&dsh ) ;
  /* Note (17Aug04):  changed '20c' to '21c' in above format code;
     The 'dsh' value was not being read correctly ... reason? (RA) */
  
  if (dsh != DSH) {  
    printf("Input grid resolution not consistent with DSH in mc_custom.h.\n Edit latter file and recompile.\n") ;
    exit(-1) ;
  }
  printf("Horizontal grid resolution: %.1f \n",dsh) ;
  for (k = 0 ; k <= kx; k++) {
    fgets(string,255,tfile) ;
    sscanf(string,"%f",&(z[k])) ;
    /* convert to km */
    z[k] /= 1000.0 ;
  }

  if (z[0] != 0.0) {
    printf("First Z-value should be zero (z[0]= %f) \n",z[0]) ;
    exit(-1) ;
  }
  x[0] = 0.0 ;
  for (i=1 ; i<=IWMAX ; i++) x[i] = x[i-1] + DSH ;
  y[0] = 0.0 ;
  for (j=1 ; j<=JWMAX ; j++) y[j] = y[j-1] + DSH ;


  /* skip four lines */
  fgets(cdum,255,tfile) ;
  fgets(cdum,255,tfile) ;
  fgets(cdum,255,tfile) ;
  fgets(cdum,255,tfile) ;

  for (k=0; k<kx; k++) {
    for (j=0; j<jx; j++) {    
      for (i=0; i<ix; i++) {
	fscanf(tfile, "%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
	       &it,&jt,&kt,
	       &(Temp[i][j][k]),
	       &absair, &abswv, &abscld,
	       &rext,&rssa,&rg, 
	       &gext,&gssa,&gg, 
	       &sext,&sssa,&sg, 
	       &iext,&issa,&ig) ;
	
	if (it != i || jt != j || kt != k) {
	  printf("Unexpected values of i,j,k encountered in input: \n(%d %d %d) vs. expected (%d %d %d)\n",
		 it,jt,kt,i,j,k) ;
	  exit(-1) ;
	}

	/* compute total extinction, SSA, and asymmetry parameter for this grid box */
	kext_total = 0.0 ;
	ksca_total = 0.0 ;
	g_total = 0.0 ;

	if (include_gas) {
	  kext_total += absair ;
	  if (include_vapor) {
	    kext_total += abswv ;
	  }
	}
	if (include_cloud) {
	  kext_total += abscld ;
	}
	if (include_rain) {
	  kext_total += rext ;
	  ksca_total += rssa*rext ;
	  g_total += rssa*rext*rg ;
	}
	if (include_graupel) {
	  kext_total += gext ;
	  ksca_total += gssa*gext ;
	  g_total += gssa*gext*gg ;
	  
	}
	if (include_ice) {
	  kext_total += iext ;
	  ksca_total += issa*iext ;
	  g_total += issa*iext*ig ;

	}
	if (include_agg) {
	  kext_total += sext ;
	  ksca_total += sssa*sext ;
	  g_total += sssa*sext*sg ;

	}
	ss = 0.0 ;
	if (kext_total > 0.0) ss = ksca_total/kext_total ;
	if (ksca_total > 0.0) 
	  g_total /= ksca_total ;
	else
	  g_total = 0.0 ;

	kext[i][j][k] = kext_total ;
	ssag[i][j][k] = ss ;
	assym[i][j][k] = g_total ;
      }
    }
  }
  fclose(tfile) ;
}


par_init() 
{
  
  printf("Enter no. of photons per pixel: ");
  scanf("%d",&nmax); 
  printf("%d\n",nmax);
  
  /* initialize vector for starting direction of photons */
  /*  theta0 = 47.0 ; */
  /*  phi0 = 35.0 ;   */
  
  /* ===============================================================
     Update (8Sep04) for purpose of Python-controlled MC model:
     - Replaced the hardcoded geometry information with
       input from the .in file; Values of theta and phi are 
       specified in the dictionary (insim.python.dict), and the 
       script assigns the proper theta for a given channel.
       The .in file (mc_directives_c.in) is updated on the fly
       by the python script 'monteCarlo'. 
       (contact Ryan_Aschbrenner@aer.com for details)
     ===============================================================
  */
    
  printf("Enter desired theta: ");
  scanf("%f",&theta0);
  printf("%f\n",theta0);
  printf("Enter desired phi: ");
  scanf("%f",&phi0);
  printf("%f\n",phi0);
  
  k10 = cos(phi0*DEG2RAD)*sin(theta0*DEG2RAD) ;
  k20 = sin(phi0*DEG2RAD)*sin(theta0*DEG2RAD) ;
  k30 = -cos(theta0*DEG2RAD) ;
  
  /* initialize all grids required by Monte Carlo model */
  init_optics() ;
  
}


write_output_file()
{
  int ipixels, jpixels, i, j ;
  float dx ;
  FILE *ofd ;


  dx = DSH/PPB ;
  ipixels = IPIXELS ;
  jpixels = JPIXELS ;

  ofd = fopen(outfile, "w") ;
  if (ofd == (FILE *) NULL) {
    printf("Couldn't open output file %s \n",outfile) ;
    exit(-1) ;
  }

  fprintf(ofd, "%.4f  =  phi0   \n", phi0 ) ;
  fprintf(ofd, "%.4f  =  theta0 \n", theta0 ) ;
  fprintf(ofd, "%.4f  =  k10    \n", k10 ) ;
  fprintf(ofd, "%.4f  =  k20    \n", k20 ) ;
  fprintf(ofd, "%.4f  =  k30    \n", k30 ) ;
  fprintf(ofd, "%.4f  =  freq   \n", freq ) ;
  fprintf(ofd, "%.4f  =  dx     \n", dx ) ;
  fprintf(ofd, "%4d  =  ipixels\n", ipixels ) ;
  fprintf(ofd, "%4d  =  jpixels\n", jpixels ) ;

  fprintf(ofd, "TBV, TBH  [0:ipixels-1][0:jpixels-1]\n") ;
  for (j=0; j<jpixels; j++) {
    for (i=0; i<ipixels; i++) {
      fprintf(ofd, "%4d %4d %5.1f %5.1f\n",i,j,tbv[i][j],tbh[i][j]) ;
    }
  } 
  
  fclose(ofd) ;
}

