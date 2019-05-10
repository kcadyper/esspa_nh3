C--------------------------------------------------------------------------

      subroutine GLIMITS(AL,FLAM,GAM,D1,D2)

C Subroutine to choose limits of integration for
C modified gamma distribution 
C FLAM is calculated assuming diameters, so this returns diameters
C D1, D2.  

      implicit none
      REAL AL,FLAM,GAM,D1,D2,VMODE,VMAX,VTARG,FACTOR,V1,V2,GAMMA_MOD
      integer IT


C     find value of diameter corresponding to volume mode
      VMODE = ((AL+3.0)/(FLAM*GAM))**(1./GAM)

      
C     find (relative) value of distribution at volume mode
      VMAX = GAMMA_MOD(AL,FLAM,GAM,VMODE)*(VMODE**3)
 
C    find lower and upper limits of integration so that kernel value
C    falls below specified threshold.

      VTARG = 1.0e-3*VMAX

C for lower limit, start at peak and then decrease by increments 
C of 20% until we get below specified threshold

      IT = 0
      FACTOR = 0.8
      D1 = VMODE
      V1 = VMAX
      do while(V1.gt.VTARG.and.IT.le.200)
         D1 = FACTOR*D1
         V1 = GAMMA_MOD(AL,FLAM,GAM,D1)*(D1**3)
         IT = IT + 1
      end do

      if(IT.ge.200)stop 'Unable to find D1.'

C for upper limit, start at peak and then inrease by increments 
C of 20% until we get below specified threshold

      IT = 0
      FACTOR = 1.2
      D2 = VMODE
      V2 = VMAX
      do while(V2.gt.VTARG.and.IT.le.200)
         D2 = FACTOR*D2
         V2 = GAMMA_MOD(AL,FLAM,GAM,D2)*(D2**3)
         IT = IT + 1
      end do

      if(IT.ge.200)stop 'Unable to find D2.'

      return
      end

C-------------------------------------------------------------------------
C Modified Gamma Distribution (sans N_0)

      function GAMMA_MOD(AL,FLAM,GAM,D)
      implicit none
      real AL,FLAM,GAM,D,GAMMA_MOD

C D is the liquid equiv. diameter of the particle [meters]
      if (D .lt. 0.0) then
         stop 'D should non-negative in Gamma_Mod()'
      else if (D .eq. 0.0 .and. AL .lt. 0.0) then
         stop 'Invalid Alpha or D in Gamma_Mod()'
      endif

      GAMMA_MOD = (D**AL)*exp(-1.0*FLAM*(D**GAM))      

      return
      end



