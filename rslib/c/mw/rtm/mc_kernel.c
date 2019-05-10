#include "mc.h"

double uran()
{
  double x ;
  static int flag=1 ;

  /* if this is the first call, set seed to known value */
  if (flag) {
    srandom(1) ;
    flag = 0 ;
  }

  /*  x = random()/2147483649.0 ;  */

  /* bias results slightly to avoid ever returning exactly 
     zero or exactly 1 */
  x = (random()+ (double) 1)/(RAND_MAX + (double) 2) ;
  return(x) ;
}

double get_rand_pathlength()
{
  return(-log(uran())) ;
}

/* start new photon at location [istart,jstart] in image array */

void init_photon(int istart,int jstart)
{
  int ipix = istart ;
  int jpix = jstart ;

  while (ipix < 0) ipix += IPIXELS ;
  while (ipix >= IPIXELS) ipix -= IPIXELS ;
  while (jpix < 0) jpix += JPIXELS ;
  while (jpix >= JPIXELS) jpix -= JPIXELS ;

  pixflag[ipix][jpix] = TRUE ;
  ncount[ipix][jpix]++ ;

/* initial direction */
  k1 = k10 ;
  k2 = k20 ;
  k3 = k30 ;

/* initial position -- allow offset so that image plane is projected onto floor of
model domain */
  d3 = z[KBMAX+1] ;
  d1 = (ipix + uran())*DPX + k1*d3/k3 ; /*DPX=(resol./pix per box)*/
  d2 = (jpix + uran())*DPY + k2*d3/k3 ; 
  if (d1 > XMAX)  /* XMAX is length of model domain in in km */
    d1 -= XMAX  ;
  else if (d1 < 0.0) 
    d1 += XMAX  ;

  if (d2 > YMAX) 
    d2 -= YMAX  ;
  else if (d2 < 0.0) 
    d2 += YMAX  ;

  /* initial grid box */
  ib = (int) (d1/DSH) ;
  jb = (int) (d2/DSH) ;
  kb = KBMAX ;
  
/* initially unscattered */
  scattered = 0 ;
  
  sigma = get_rand_pathlength() ;   /* optical distance to next scattering */
  tauh = 1.0 ;                    /* attenuation experienced thus far */
  tauv = 1.0 ;
  tbphoth = 0.0 ;                 /* cumulative temperature of this photon */
  tbphotv = 0.0 ;
}

void get_scattering_pars();
scatter()
{
  float deltb, em ;

  scattered++ ;
  get_scattering_pars() ;   /* get new direction and single-scatter albedo */
  em = 1.0 - ssa ;
  deltb = em*local_temp() ;
  tbphoth += deltb*tauh ;
  tbphotv += deltb*tauv ;
  tauh = tauh*ssa ;
  tauv = tauv*ssa ;
  sigma = get_rand_pathlength() ;   /* optical distance to next scattering */
}

find_box_exit()
{
  int id,jd,kd ;
  int iw,jw,kw ;
  float s1,s2,s3 ;
  
  if (k1 == 0.0) id = -1 ;
  else { if (k1 > 0.0) id = 1; else id = 0 ; }
  if (k2 == 0.0) jd = -1 ;
  else { if (k2 > 0.0) jd = 1; else jd = 0 ; }
  if (k3 == 0.0) kd = -1 ;
  else { if (k3 > 0.0) kd = 1; else kd = 0 ; }

  if (id >= 0) {       /* determine distance to each coordinate plane */
    iw = ib + id ;
    s1 = (x[iw] - d1)/k1 ;
  } else s1 = BIG ;   /* "infinity" */
  
  if (jd >= 0) {
    jw = jb + jd ;
    s2 = (y[jw] - d2)/k2 ;
  } else s2 = BIG ;
  
  if (kd >= 0) {
    kw = kb + kd ;
    s3 = (z[kw] - d3)/k3 ;
  } else s3 = BIG ;
  
  if (s1 < s2) {
    if (s1 < s3) {
      s = s1 ;
      d1new = x[iw] ;
      d2new = d2 + s*k2 ;
      d3new = d3 + s*k3 ;
      if (id == 0) ibnew = ib - 1 ; else ibnew = ib + 1 ;
      jbnew = jb ;
      kbnew = kb ;
    } else {
      s = s3 ;
      d1new = d1 + s*k1 ;
      d2new = d2 + s*k2 ;
      d3new = z[kw] ;
      ibnew = ib ;
      jbnew = jb ;
      if (kd == 0) kbnew = kb - 1 ; else kbnew = kb + 1 ;
    }
  } else {
    if (s2 < s3) {
      s = s2 ;
      d1new = d1 + s*k1 ;
      d2new = y[jw] ;
      d3new = d3 + s*k3 ;
      ibnew = ib ;
      if (jd == 0) jbnew = jb - 1 ; else jbnew = jb + 1 ;
      kbnew = kb ;
    } else {
      s = s3 ;
      d1new = d1 + s*k1 ;
      d2new = d2 + s*k2 ;
      d3new = z[kw] ;
      ibnew = ib ;
      jbnew = jb ;
      if (kd == 0) kbnew = kb - 1 ; else kbnew = kb + 1 ;
    }
  }
}

void get_reflection_pars();
lower_boundary()
{

  float em, tsurf ;

  kb = 0 ;
  get_reflection_pars() ;        /* get new direction and reflectivity */
  tsurf = surface_temp() ;
  em = 1.0 - reflecth ;
  tbphoth += tauh*em*tsurf ;
  em = 1.0 - reflectv ;
  tbphotv += tauv*em*tsurf ;
  tauh = tauh*reflecth ;
  tauv = tauv*reflectv ;
  sigma = get_rand_pathlength() ;   /* optical distance to next scattering */
}

upper_boundary()
{
  float c ;

  tbphotv += tauv*TCOSMIC ;         /* assumes that top level corresponds */
  tbphoth += tauh*TCOSMIC ;         /*     to top of atmosphere (temporary)*/

  tauv = 0.0 ;   
  tauh = 0.0 ;
}

move_photon()
{
  float dsigma ;

/* get coordinates of exit point (d1new,d2new,d3new) and
   of next box  (ibnew,jbnew,kbnew) */
  find_box_exit() ;  

/* optical distance to exit */
  dsigma = s*kext[ib][jb][kb] ;

/* if sigma large enough, then photon reaches exit point */
  if (dsigma < sigma) {
    ib = ibnew; jb = jbnew; kb=kbnew;
    d1 = d1new; d2 = d2new; d3=d3new;
    sigma = sigma - dsigma ;

/* here we assume periodic lateral boundary conditions for domain */
    if (ib < 0) {
      ib = IBMAX ;
      d1 = XMAX ;
    } else if (ib > IBMAX) {
      ib = 0 ;
      d1 = 0.0 ;
    }
    if (jb < 0) {
      jb = JBMAX ;
      d2 = YMAX ;
    } else if (jb > JBMAX) {
      jb = 0 ;
      d2 = 0.0 ;
    }

/* check for photon encountering upper or lower boundary of domain */
    if (kb < 0) 
      lower_boundary() ;
    else if (kb > KBMAX)
      upper_boundary() ;

/* otherwise photon stays in current box and experiences extinction event */
  } else {
    s = sigma/kext[ib][jb][kb] ;
    d1 = d1 + s*k1 ;
    d2 = d2 + s*k2 ;
    d3 = d3 + s*k3 ;
    scatter() ;
  }
}

tally_photon(int istart, int jstart)
{
  int ipix = istart ;
  int jpix = jstart ;

  while (ipix < 0) ipix += IPIXELS ;
  while (ipix >= IPIXELS) ipix -= IPIXELS ;
  while (jpix < 0) jpix += JPIXELS ;
  while (jpix >= JPIXELS) jpix -= JPIXELS ;

  if (tauh < 1.0) tbphoth = tbphoth/(1.0-tauh) ;
  if (tauv < 1.0) tbphotv = tbphotv/(1.0-tauv) ;
  tbh[ipix][jpix] += tbphoth ;
  tbv[ipix][jpix] += tbphotv ;
}

void write_output_file();
final_tally() 
{ 
  float tv,th ;
  int ipixcount,nn ;
  int i,j ;

  tbhave = 0.0 ;
  tbvave = 0.0 ;
  ipixcount = 0 ;

  for (i=0; i<IPIXELS; i++) { 
    for (j=0; j<JPIXELS; j++) { 
      if (pixflag[i][j]) {
	nn = ncount[i][j] ;
	tbh[i][j] /= nn ;
	tbv[i][j] /= nn ;

	tv = tbv[i][j] ;
	th = tbh[i][j] ;

	tbhave += th ;
	tbvave += tv ;

	ipixcount++ ;
	
      } else {
	tv = 0.0 ;
	th = 0.0 ;
	tbv[i][j] = tv ;
	tbh[i][j] = th ;
      }
    }
  }
  tbhave /= ipixcount ;
  tbvave /= ipixcount ;

  write_output_file() ;
}


grid_init()
{
  int i,j ;

  for (i=0; i<IPIXELS; i++) { 
    for (j=0; j<JPIXELS; j++)	{ 
      tbh[i][j] = 0.0 ;
      tbv[i][j] = 0.0 ;
      pixflag[i][j] = FALSE ;
      ncount[i][j] = 0 ;
    } 
  } 
}


main_loop(int nmax) 
{

  int istart, jstart, n ,idum ;
  float clock_sum, old_clock, new_clock, diff_clock ;
  
/* zero out results grids */
  grid_init() ;

  clock_sum = 0.0 ;
  old_clock = ((float) clock())/CLOCKS_PER_SEC ;

  printf("***BEGIN MAIN LOOP\n") ;

  /* release NMAX photons for each pixel */
  printf("%2i IPIXELS: ",IPIXELS,JPIXELS);
  for (n=0 ; n < nmax ; n++) { 
    if ((n % 50) == 0) {  
      /* 'if' statement added to inform of model progress (RA) */
      printf("%.3f percent complete\n",n/(float)nmax);
    }
    /*printf("%i .. \n",n);*/
    /* loop over image pixels */
    for (istart = 0 ; istart < IPIXELS ; istart++) {
      for (jstart = 0 ; jstart < JPIXELS ; jstart++) { 
	if ((istart == 0 && jstart == 7) || (istart == 53 && jstart == 5)) {
	  idum = 1 ;
	  idum = 2 ;
	}
	init_photon(istart,jstart) ;
	while ((tauv > MINTAU)||(tauh > MINTAU)) {
	  move_photon() ;
	}
	tally_photon(istart,jstart) ;
	/*if (istart == 0 && jstart < 3) {
	  printf("%.3f %.3f \n",tbh[ipix][jpix],tbv[ipix][jpix]) ;
	  }*/
      }
    }
    new_clock = ((float) clock())/CLOCKS_PER_SEC ;
    diff_clock = new_clock - old_clock ;
    old_clock = new_clock ;
    if (diff_clock > 0.0) clock_sum += diff_clock ;
  }
  printf("***END MAIN LOOP: ELAPSED CPU TIME = %f\n", clock_sum) ;

/* process results */
  final_tally() ;

/* print domain-averaged results to console */

  printf("TV = %5.1f TH = %5.1f \n",tbvave,tbhave) ;
  
}



