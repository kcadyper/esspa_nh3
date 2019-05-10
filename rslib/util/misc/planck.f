c-------------------------------------------------------------------------------
c     $Name$ 
c     $Id$ 
c     Copyright AER, Inc., 2002, 2003. All rights Reserved.
c-------------------------------------------------------------------------------

      function bt(vn,rad)
c*******************************************************************
c* Function name: bt
c* Purpose: Converts radiance to brightness temperature
c* Usage: brightTemp = bt(vn, rad)
c* Description: At specified frequency, the function converts
c*              the radiance (mw/m2/str/cm-1) to the brightness 
c*              temperature (K).
c* Inputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* vn           real    frequency (wave number)
c* rad          real    radiance (mw/m2/str/cm-1)
c* 
c* Outputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* bt           real    brightness temperature (K)
c* 
c* Common blocks: none
c* Includes: none
c* Externals: none
c*
c* Copyright: Atmospheric and Environmental Research, Inc., 1997        
c* Developed by Atmospheric and Environmental Research, Inc.            
c*******************************************************************
c $ radiance to brightness temperature (hmw)
c * 'new' planck's constant, velocity of light, boltzmann's constant

      real rad,bt,vn
      parameter (h = 6.626176e-27, c = 2.997925e+10, b = 1.380662e-16)
      parameter (c1 = 2.e0*h*c*c)
      parameter (c2 = h*c/b)

      tnb(x,y,z)=y/alog(x/z+1.e0)
c
      f1=c1*vn**3
      f2=c2*vn
      r=rad
      tbb=tnb(f1,f2,r)
      bt=tbb
      return
      end

      function dbdt(vn,t)
c*******************************************************************
c* Function name: dbdt
c* Purpose: Calculates the derivative of radiance w.r.t. temperature
c* Usage: dRaddTemp = dbdt(vn, t)
c* Description: At specified frequency, the function calcualtes
c*              the derivative of the radiance (or brightness 
c*              temperature) w.r.t. temperature (K). This is a 
c*              double-precision function. Please refer to draddt 
c*              for single-precision process.
c* Inputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* vn           real*8  frequency (wave number)
c* t            real*8  temperature (K)
c* 
c* Outputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* dbdt         real*8  derivative of brightness temperature
c*                      w.r.t. temperature
c* 
c* Common blocks: none
c* Includes: none
c* Externals: none
c*
c* Copyright: Atmospheric and Environmental Research, Inc., 1997        
c* Developed by Atmospheric and Environmental Research, Inc.            
c*******************************************************************
      implicit real*8 (a-h,o-z)
      real t,dbdt,vn
      parameter(h=6.626176d-27,c=2.997925d+10,b=1.380662d-16)
      parameter(c1=2.d0*h*c*c)
      parameter(c2=h*c/b)
c
      t_1=1/t
c
      ct1v1=c1*vn**3
      ct2v1=c2*vn*t_1
      ct3v1=exp(ct2v1)
c
      bplan=ct1v1/(ct3v1-1.)
      dbdt=bplan/(ct3v1-1.)*ct2v1*ct3v1*t_1
c
      return
      end

      function draddt(vn,t)
c*******************************************************************
c* Function name: draddt
c* Purpose: Calculates the derivative of radiance w.r.t. temperature
c* Usage: dRaddTemp = draddt(vn, t)
c* Description: At specified frequency, the function calcualtes
c*              the derivative of the radiance (or brightness 
c*              temperature) w.r.t. temperature (K). This is 
c*              a single-precision function. Please refer to 
c*              dbdt for dboule-precision process.
c* Inputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* vn           real    frequency (wave number)
c* t            real    temperature (K)
c* 
c* Outputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* draddt       real    derivative of radiance (or brightness
c*                      temperature) w.r.t. temperature
c* 
c* Common blocks: none
c* Includes: none
c* Externals: none
c*
c* Copyright: Atmospheric and Environmental Research, Inc., 1997        
c* Developed by Atmospheric and Environmental Research, Inc.            
c*******************************************************************
      real t,draddt,vn
      parameter(h=6.626176e-27,c=2.997925e+10,b=1.380662e-16)
      parameter(c1=2.e0*h*c*c)
      parameter(c2=h*c/b)
c
      t_1=1/t
c
      ct1v1=c1*vn**3
      ct2v1=c2*vn*t_1
      ct3v1=exp(ct2v1)
c
      bplan=ct1v1/(ct3v1-1.)
      draddt=bplan/(ct3v1-1.)*ct2v1*ct3v1*t_1
c
      return
      end
                                                                   

      function wnplan(vn,tem)
c*******************************************************************
c* Function name: wnplan
c* Purpose: Double-precision planck function
c* Usage: wn = wnplan(vn, tem)
c* Description: At specified frequency, the function calcualtes
c*              the radiance at given temperature using the
c*              planck equation. Please refer to planck for
c*              the single-precision process.
c* Inputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* vn           real*8  frequency (wave number)
c* tem          real*8  temperature (K)
c* 
c* Outputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* wnplan       real*8  radiance (mw/m2/str/cm-1)
c* 
c* Common blocks: none
c* Includes: none
c* Externals: none
c*
c* Copyright: Atmospheric and Environmental Research, Inc., 1997        
c* Developed by Atmospheric and Environmental Research, Inc.            
c*******************************************************************
c $ temperature to planck radiance (hmw)
c $ 'new' planck's constant, velocity of light, boltzmann's constant
c units are mw*cm/m2/sr
      implicit real*8 (a-h,o-z)
      real tem,wnplan,vn
      parameter (h = 6.626176d-27, c = 2.997925d+10, b = 1.380662d-16)
c units of h are mw, units of c are cm/s
      parameter (c1 = 2.d0*h*c*c)
      parameter (c2 = h*c/b)
      bnt(x,y,z)=x/(exp(y/z)-1.d0)
c
      f1=c1*vn**3
      f2=c2*vn
      t=tem
      rad=bnt(f1,f2,t)
      wnplan=rad
      return
      end

      function planck(vn,tem)
c*******************************************************************
c* Function name: planck
c* Purpose: Single-precision planck function
c* Usage: plnk = planck(vn, tem)
c* Description: At specified frequency, the function calcualtes
c*              the radiance at given temperature using the
c*              planck equation. Please refer to wnplan for
c*              the double-precision process.
c* Inputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* vn           real    frequency (wave number)
c* tem          real    temperature (K)
c* 
c* Outputs:
c* Var_name     Type    Description
c* --------     ----    -----------
c* planck       real    radiance (mw/m2/str/cm-1)
c* 
c* Common blocks: none
c* Includes: none
c* Externals: none
c*
c* Copyright: Atmospheric and Environmental Research, Inc., 1997        
c* Developed by Atmospheric and Environmental Research, Inc.            
c*******************************************************************
c $ temperature to planck radiance (hmw)
c $ 'new' planck's constant, velocity of light, boltzmann's constant
c units are mw*cm/m2/sr
      real tem,planck,vn
      parameter (h = 6.626176e-27, c = 2.997925e+10, b = 1.380662e-16)
c units of h are mw, units of c are cm/s
      parameter (c1 = 2.e0*h*c*c)
      parameter (c2 = h*c/b)
      bnt(x,y,z)=x/(exp(y/z)-1.e0)
c
      f1=c1*vn**3
      f2=c2*vn
      t=tem
      rad=bnt(f1,f2,t)
      planck=rad
      return
      end

