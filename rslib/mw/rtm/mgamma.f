      subroutine MGAMMA(F1,F2,TC,NU,FLAG,AL,FLAM,GAM,MGMEXT
     &,MGSSA,MGASYMM,D1,D2)

C variable declarations
      implicit none
      integer NANG, FLAG
      real MGMEXT,ASYMM,X,x2,x4
      real EXTCS,MGASYMM,VOL,DAREA,DVOL,SCACS,PI,D,CC
      real QEXT,QSCA,F1,F2,NU,TC,RHOAV,K2,GSCA,QBK
      real COUNT,MGSSA,RATIO,RHOW,RHOI,DPHY,WAVEL
      real GAMMA_MOD,AL,FLAM,GAM,D1,D2,DSTP,G,S1(1999),S2(1999)
      complex NAV,m,ac,bc,m2,c1,c2

C*************************************************************************
C
C     Subroutine MGAMMA.F provides the calling program or routine with 
C     the single scatter albedo, asymmetry parameter, and mass extinction
C     coefficient for a distribution of spherical hydrometeors.  Each 
C     diameter is treated, and the results are summed within the main loop.
C     The distribution is a modified Gamma distribution, see GLIMITS.F 
C     GAMMA_MOD fuction for the form of the distribution.
C     
C     INPUT : The subroutine accepts (and requires) the following inputs
C             F1 = Fraction (out of 1) of liquid water.
C             F2 = Fraction (out of 1, but not greater than 1-F1) of ice.
C             TC = Temperature of the system [celcius]
C             NU = Frequency of incident radiation [gigahertz]
C             FLAG = value to choose which mixing formula to use
C                    using Bruggeman (MIXFLAG = 0)
C                    MIXFLAG for Maxwell Garnet would be 
C                    1230  for [[[Water(1)], Ice(2)], Air(3)]
C                    1231  for [Water(1),[[Ice(2)], Air(3)]]]
C                    1320  for [[[Water(1)], Air(3)], Ice(2)]
C 
C             AL = alpha parameter (see GLIMITS.F GAMMA_MOD function)
C             FLAM = lambda parameter (see GLIMITS.F GAMMA_MOD function)
C             GAM = gamma parameter (see GLIMITS.F GAMMA_MOD function)
C             D1 = Minimum liquid equivalent diameter (obtained from 
C                  GLIMITS.F) [meters]
C             D2 = Maximum liquid equivalent diameter (obtained from 
C                  GLIMITS.F) [meters]
C     
C     OUTPUT : The subroutine returns the following values to the calling
C             MGMEXT = Mass extinction (attenuation) coefficient [m^2/kg]
C             MGSSA = Single scatter albedo (as defined in Bohren Huffman)
C             MGASYMM = Asymmetry parameter (cosine average of scattering 
C                     angle)  g = <cos(theta)>
C
C     NOTES  :
C             This subroutine is most effective when compiled
C             using the +autodblpad flag with F77.  This
C             allows for values of alpha higher than ~75.
C             
C             This subroutine calls the following routines/functions:
C             RHOW(TC) = density of water at temperature TC [gm/cm^3]
C             RHOI(TC) = density of ice at temperature TC   [gm/cm^3]
C             GAMMA_MOD(AL,FLAM,GAM,D) = value of the modified gamma 
C               distribution for a specific diameter D.
C
C     Other Variables:
C             VOL,EXTCS,SCACS,ASYMM are summation variables
C             G = Value of the mod. gamma distribution for a given value
C                 of diameter D
*
* revision history
* 12/28/2000   added small-particle approximation (G. Petty)
C*************************************************************************


C parameters
      PI = 4.E0*ATAN(1.E0)
      CC = 2.99792E+8
      NANG = 1

C initialize variables
      VOL = 0.
      D = 0.
      COUNT = 0.
      EXTCS = 0.
      SCACS = 0.
      ASYMM = 0.

C Calculate the average density

      RHOAV = (F1*RHOW(TC) + F2*RHOI(TC))*1.0e3

C calculate ratio of diameter to standardized liquid equivalent diameter at 0 C

      RATIO = (RHOW(0.0)*1.0e3/RHOAV)**(1./3.)

C Get discretization interval
      DSTP = ((D2 - D1)/1000.0)
      if (DSTP.eq.0) pause 'DSTP = 0'

C Begin integrating
      D = D1

C open a temporary file to write gamma dist data to
C      open(unit=2, file='temp.out', status = 'unknown')

C ------- begin main loop --------
      do while(D.le.D2)

C Diameter : D [meters] 
C Convert to Physical Diameter from Liquid Equivalent Diameter         
C DPHY is in meters, D is in meters
        DPHY = (D*RATIO)

C Area : DIAMETER squared [m^2]

        DAREA = (PI/4.)*DPHY**2


C Volume : DIAMETER cubed [m^3]

        DVOL = (PI/6.)*DPHY**3


C Get Gamma distribution value given the liquid equivalent
C diameter D [meters]

        G = GAMMA_MOD(AL,FLAM,GAM,D)

C NU MUST BE IN GHZ WHEN PASSED TO CDMX

        call CDMX(F1,F2,TC,NU,FLAG,K2,NAV)
        WAVEL= (CC/(NU*1.E9))
C Size parameter for diameter (factor of 2 smaller than eqn for radius)
        X    = PI*DPHY/WAVEL

* use small particle approximation for small x

        if (x .lt. 0.3) then
           m = nav
           m2 = m**2
           x2 = x*x
           x4 = x2*x2

* from Bohren and Huffman, Eqs. 5.7-5.9

           ac = (m2 - 1)/(m2 + 2)
           bc = (m2**2 + 27.*m2 + 38)/(2*m2 + 3)
           c1 = 4*x*aimag(ac*(1. + (x2/15.)*ac*bc))
           c2 = (8./3.)*x4*real(ac*ac)
           qext = c1 + c2
           qsca = (8./3.)*x4*real(cabs(ac)**2)
           qbk = qsca*(3./8.)*4.
           gsca = 0.0

* delete following line after testing
*           call BHMIE(X,NAV,NANG,S1,S2,QEXT,QSCA,QBK,GSCA)

        else
           call BHMIE(X,NAV,NANG,S1,S2,QEXT,QSCA,QBK,GSCA)
        endif

C sum up various radiative cross-sections
        EXTCS = EXTCS + QEXT*DAREA*G
        SCACS = SCACS + QSCA*DAREA*G
        ASYMM = ASYMM + QSCA*DAREA*GSCA*G

C sum up volume, increment diameter etc
        VOL = VOL + DVOL*G
        D = D + DSTP
C        COUNT = COUNT + 1

C have we finished integrating?  If not, continue loop

      end do

C we're out of the loop; now consolidate results
C mass extinction coeff, single scatter albedo, asymmetry param

      MGMEXT  = EXTCS/(VOL*RHOAV)
      MGSSA   = SCACS/EXTCS
      MGASYMM = ASYMM/SCACS

      return
      end





