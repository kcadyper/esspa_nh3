#define HENYEYGREENSTEIN
#undef ISOTROPIC 

#undef LAMBERTIAN
#define SPECULAR


/* gridbox definitions */

#define IBOXES 66
#define JBOXES 66
#define KBOXES 30

#define KPARLEVELS (KBOXES + 1)

/* horizontal resolution (km) */
#define DSH 3.0

/* dimensions of domain (km) */
#define XMAX (IBOXES*DSH)
#define YMAX (JBOXES*DSH)

#define PPB 1               /* image pixels per horizontal box */






































