      subroutine MC_GAMMA(PFLAG,W,F1,F2,DD,PR,TK,NU,MGSSA,MGMEXT
     &,MGASYM)

C Heavily modified version of plotmg for use with monte_carlo code.
C 10/14/1999 BTJ

C Inputs : F1,F2,TC,NU
C Outputs: MGSSA,MGMEXT,MGASYM

      implicit none
      real F1,F2,TC,NU,A1,A2,W,PR,TK,DD,RHOI
      real AL,FLAM,GAM,MGSSA,MGMEXT,MGASYM,RHOG
      integer MIXFLAG,PFLAG

C set modified gamma distribution paramters such that we have an
C exponential distribution.

      AL=0.0
      GAM=1.0

      if (PFLAG.eq.0) then
          call FLAM_RAIN(W,DD,PR,TK,FLAM)
      end if
      if (PFLAG.eq.1) then
         call FLAM_SNOW(W,DD,PR,TK,FLAM)
      end if
      if (PFLAG.eq.2) then
         TC = TK - 273.15
         RHOG = 1000*F2*RHOI(TC)
         call FLAM_GRAUPEL(W,DD,RHOG,PR,TK,FLAM)
      end if   
      write(*,*)"FLAM:", FLAM
C get limits of integration A1, A2 [meters]

      call GLIMITS(AL,FLAM,GAM,A1,A2)

C N(D)=D^AL*exp(-FLAM*D^GAM)

      TC = TK-273.15

C using Bruggeman (MIXFLAG = 0)
C MIXFLAG for Maxwell Garnet would be 
C 1230  for [[[Water(1)], Ice(2)], Air(3)]
C 1231  for [Water(1),[[Ice(2)], Air(3)]]]
C 1320  for [[[Water(1)], Air(3)], Ice(2)]
C etc. 

      MIXFLAG = 0

C Inputs F1,F2,TC,NU,MIXFLAG,AL,FLAM,GAM,A1,A2
C Outputs MGMEXT [m^2/kg],MGSSA,MGASYM
      call MGAMMA(F1,F2,TC,NU,MIXFLAG,AL,FLAM,GAM,MGMEXT
     & ,MGSSA,MGASYM,A1,A2)

      return
      END








