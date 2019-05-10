c &"$Id$"/

      SUBROUTINE FIND_EMISS_TOT (FREQ,THT,SST,WIND,SAL,E1,E2)
c
c     NAME: 
c     FIND_EMISS_TOT
c     
c     USAGE:
c     CALL FIND_EMISS_TOT (FREQ,THT,SST,WIND,EMISS) 
c
c     DESCRIPTION:
c     calculates isotropic part of total sea surface emissivity for V and H pol 
c     for given FREQ, EIA, SST and wind speed

c 
c
c
c
c    INPUT:
c    NAME       DESCRIPTION                 TYPE                   LENGTH      UNIT       RANGE
c ----------------------------------------------------------------------------------------------
c
c    FREQ       FREQUENCY                          REAL           SCALAR  GHz       > 0
c    THT        EIA at IFREQ               REAL           SCALAR  deg       [ 0,88]
c    SST        Sea Surface Temperature    REAL           SCALAR  Celsius   [-2,40]
c    WIND       Ocean Surface Wind Speed   REAL           SCALAR  m/s       [ 0,40]
c                   10 m above ground                 
c    SAL        SALINITY                   REAL           SCALAR  ppt       [ 0,40]           
c
c    OUTPUT:
c    NAME       DESCRIPTION                           TYPE           LENGTH           RANGE
c ----------------------------------------------------------------------------------------------
c    EMISS     total sea surface emissivity          VECTOR                  (2)                [0,1]  
c    E(W,T) (1:V  2:H) 
c



      IMPLICIT NONE
      REAL :: FREQ,THT,SST,WIND,SAL,e1,e2
      REAL, DIMENSION(2) :: E0,DE,EMISS
      SAVE

      CALL  FINDEM0(FREQ,THT,SST,SAL, E0)   ! specular
      CALL  FIND_WIND_SIGNAL (FREQ,THT,SST,WIND,  DE)  ! wind induced
      EMISS = E0 + DE  ! total
        e1=EMISS(1)
        e2=EMISS(2)
      RETURN
      END


      subroutine FINDEM0(FREQ,THT,SST,SAL, E0)
c     specular emissivity for v and h-pol

C   INPUT:
c   NAME   PARAMETER  UNIT  RANGE
c   FREQ   FREQUENCY  [GHz] 1 to 400 (?)
c   SST    SST        [Celsius   -25 C  to 40 C for pure water
c                                 -2 C  to 40 C for saline water
c   SAL    SALINITY   [ppt]  0 to 40
c
c   OUTPUT:
c   E0(1:2)     SPECULAR SEA SURFACE EMISSIVITY (1= v-pol,  2= h-pol)


      implicit none
      REAL(4) :: FREQ,THT,SST,SURTEP,SAL,E0(2)
      COMPLEX(4) :: PERMIT

      SURTEP = SST + 273.15
      CALL meissner_wentz(FREQ,SURTEP,SAL, PERMIT)
      CALL FRESNEL(PERMIT,THT, E0)

      RETURN
      END



c    DIELECTRIC CONSTANT
      SUBROUTINE  meissner_wentz(freq,t,s,eps)
C    COMPLEX DIELECTRIC CONSTANT: EPS
c    T. Meissner, February 2002

C   INPUT:
c   NAME   PARAMETER  UNIT  RANGE
c   FREQ   FREQUENCY  [GHz] 1 to 400
c   T      SST        [K]   248.15 K (-25 C) to 313.15 K (40 C) for pure water
c                           271.15 K (-2  C) to 313.15 K (40 C) for saline water
c   S      SALINITY   [ppt]  0 to 40
c
c   OUTPUT:
c   EPS    COMPLEX DIELECTRIC CONSTANT 
c          NEGATIVE IMAGINARY PART TO BE CONSISTENT WITH WENTZ1 CONVENTION
c
c
c
      implicit none
    
      REAL(4) :: FREQ,T,SST,S
      REAL(4) :: E0,E1,E2,N1,N2
      REAL(4) :: A0,A1,A2,B1,B2
      REAL(4) :: E0S,E1S,E2S,N1S,N2S
      REAL(4) :: SIG35,R15,RTR15,ALPHA0,ALPHA1,SIG,F0=17.97510
      COMPLEX(4) :: J = (0.0,1.0), EPS
      INTEGER :: ISTART = 0 
c     EXTERNAL ASCII FILE 
c     PATH + FILENAME NEED TO BE SPECIFIED BY USER    
      character (LEN = 100) :: inputfile
      REAL, DIMENSION(11) ::  X
      REAL, DIMENSION(13) ::  Z
      SAVE X,Z

      inputfile=
     &'/npoess/sap/sw_share/mw_ocean_model/rss_emiss/rss_dielectric.lis'
      If (ISTART == 0) then
      open(unit=3,file=inputfile,status='old',form='formatted')
      READ(3,*) X
      READ(3,*) Z(1:11)
      READ(3,*) Z(12:13)
      CLOSE(3)
      endif
      SST = T - 273.15 ! [Celsius]


c     PURE WATER

      E0    = (3.70886E4 - 8.2168E1*SST)/(4.21854E2 + SST) ! Stogryn et al.
      E1    = X(1) + X(2)*SST + X(3)*SST**2
      N1    = (45.00 + SST)/(X(4) + X(5)*SST + X(6)*SST**2)
      E2    = X(7) + X(8)*SST
      N2    = (45.00 + SST)/(X(9) + X(10)*SST + X(11)*SST**2)



c     Saline Water
c     Conductivity [S/m] taken from Stogryn et al. 
      SIG35 = 2.903602 + 8.60700E-2*SST + 4.738817E-4*SST**2 - 
     &      2.9910E-6*SST**3 + 4.3047E-9*SST**4
      R15   = S*(37.5109+5.45216*S+1.4409E-2*S**2)/
     &      (1004.75+182.283*S+S**2)

      alpha0 = (6.9431+3.2841*S-9.9486E-2*S**2)/(84.850+69.024*S+S**2)
      alpha1 = 49.843 - 0.2276*S + 0.198E-2*S**2
      RTR15 = 1.0 + (SST-15.0)*ALPHA0/(ALPHA1+SST)

      SIG = SIG35*R15*RTR15



c     Permittivity
      A0  = exp(Z(1)*S + Z(2)*S**2 + Z(3)*S*SST)
      E0S = A0*E0
      B1  = 1.0 + S*(Z(4) + Z(5)*SST + Z(6)*SST**2)
      N1S = N1*B1
      A1  = exp(Z(7)*S + Z(8)*S**2 + Z(9)*S*SST)
      E1S = E1*A1
      B2 = 1.0 + S*(Z(10) + Z(11)*SST)
      N2S = N2*B2
      A2 = 1.0  + S*(Z(12) + Z(13)*SST)
      E2S = E2*A2


c     Debye Law (2 relaxation wavelengths)
      EPS = (E0S - E1S)/(1.0 - J*(FREQ/N1S)) + 
     1      (E1S - E2S)/(1.0 - J*(FREQ/N2S)) + E2S + 
     2      J*SIG*F0/FREQ


      EPS = CONJG(EPS)

      ISTART = 1


      RETURN 
      END




      subroutine  fresnel(PERMIT,THT,E0)
c    calculates specular emissivity from complex dielectric constant
c    using FRESNEL formulas

      IMPLICIT NONE
        REAL :: DEG2RAD 
      REAL :: THT, COSTHT,SINSQTHT
      COMPLEX(4) :: PERMIT,RV,RH,ESQRT
      REAL, DIMENSION(2) :: E0
      SAVE
      DEG2RAD=acos(-1.0)/180.0
      COSTHT=COS(THT*DEG2RAD)
      SINSQTHT=1.-COSTHT*COSTHT

      ESQRT=CSQRT(PERMIT-SINSQTHT)
      RH=(COSTHT-ESQRT)/(COSTHT+ESQRT) 
      RV=(PERMIT*COSTHT-ESQRT)/(PERMIT*COSTHT+ESQRT)

      E0(1)  =1.-RV*CONJG(RV) 
      E0(2)  =1.-RH*CONJG(RH)

      RETURN
      END





      SUBROUTINE FIND_WIND_SIGNAL (XFREQ,XTHT,XSST,XWIND,DEMISS)
c     VERSION 3: JAN 2003: assume logarithmic frequency dependence
c     V-pol Zero at 85.5 GHz
c
c     NAME: 
c     FIND_WIND_SIGNAL
c     
c     USAGE:
c     CALL FIND_WIND_SIGNAL(XFREQ,XTHT,XSST,XWIND,   DEMSISS) 
c
c     DESCRIPTION:
c     calculates isotropic part of wind sea surface emissivity signal for V and H pol 
c     for given EIA, SST, FREQ, wind speed
c     logarithmic frequqnecy dependence (DECEMBER 2002)

c 
c
c
c
c    INPUT:
c    NAME       DESCRIPTION                 TYPE                   LENGTH      UNIT       RANGE
c ----------------------------------------------------------------------------------------------
c
c    XFREQ       FREQUENCY                      REAL           SCALAR  GHz          > 0
c    XTHT        EIA at IFREQ               REAL           SCALAR  deg      [0,88]
c    XSST        Sea Surface Temperature    REAL           SCALAR  Celsius    [-2,40]
c    XWIND       Ocean Surface Wind Speed   REAL           SCALAR  m/s        [0,40]
c                    10 m above ground                 
c
c    OUTPUT:
c    NAME       DESCRIPTION                           TYPE           LENGTH           RANGE
c ----------------------------------------------------------------------------------------------
c    DEMISS     isotropic wind induced emissisivities REAL           VECTOR(2)            [0,1]  
c    = E(W,T) - E0(W,T) for V and H  (1:V  2:H) 
c
 


      IMPLICIT NONE

      REAL, PARAMETER, DIMENSION(5) :: ZFREQ = 
     &     (/6.9250,10.650,18.70,36.50,85.500/)  
c     1    (/6.9250,10.650,18.70,36.50,85.500/)  
      ! frequency interpolation grid points 



      INTEGER :: IFREQ
      REAL(4) :: DF,DF0
        INTEGER, PARAMETER :: NTHT = 45,  NW = 41 , NSST = 43
      REAL, PARAMETER :: DWIND = 1.0, DSST = 1.0, DTHT = 2.0 
      REAL(4), PARAMETER :: THT0 = 0.0, SST0 = -2.0, WIND0 = 0.0


      !REAL(4), DIMENSION(5,2,NTHT,NSST,NW), SAVE :: EMISS_ARRAY
      REAL(4), DIMENSION(5,2,NTHT,NSST,NW) :: EMISS_ARRAY


      INTEGER :: IPOL,ITHT,ISST,IWIND,L
      INTEGER :: JFREQ,JPOL,JTHT,JSST

      REAL(4) :: THT,SST,WIND,FREQ
      REAL(4) :: XFREQ,XTHT,XSST,XWIND

      REAL(4), DIMENSION(2) :: DEMISS
      REAL(4) :: T_THT, T_SST, T_WIND, THT_1, SST_1, WIND_1


      INTEGER :: i1,i2,j1,j2,k1,k2
      REAl(4), DIMENSION(2,2) :: FARAY
      REAl(4), DIMENSION(2  ) :: GARAY

      REAL(4), DIMENSION(2,0:1) :: YY


      INTEGER, SAVE :: ISTART = 0!,IFREQ
      REAL(4), SAVE :: FSAVE = 0.0!,DF,DF0
      !SAVE

c     EXTERNAL ASCII FILE 
c     PATH + FILENAME NEED TO BE SPECIFIED BY USER     
      character (LEN=100) :: inputfile= 
     &   '/npoess/sap/sw_share/mw_ocean_model/rss_emiss/emiss_wind.lis'


c     WHEN THE SUBROUTINE IS CALLED THE FIST TIME (ISTART = 0), EMISS_ARRAY 
c     are read from inputfile
c     for subsequent calls istart = 1 and the arrays from the last run are used  

      if (istart == 0) then
      open(unit=3,file=inputfile,form='formatted',status='old')

      DO JFREQ = 1,1
      DO JPOL  = 1,2
      DO JTHT  = 1,NTHT
      DO JSST  = 1,NSST
        !print *,JFREQ,JPOL,JTHT,JSST
      READ(3,*) EMISS_ARRAY(JFREQ,JPOL,JTHT,JSST,1:NW)
        !print *,nw,EMISS_ARRAY(JFREQ,JPOL,JTHT,JSST,1:NW)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      close(3)
      istart = 1
      endif


      freq = xfreq

c    frequency interpolation points
      if (abs(FREQ-FSAVE) > 1.0E-3) then
      if (FREQ < ZFREQ(1)) then
         IFREQ = 1
      ELSE if (FREQ >= ZFREQ(5)) then
         IFREQ = 4
      ELSE  
      do JFREQ = 1,4
      if (freq >= zfreq(jfreq) .and. freq < zfreq(jfreq+1) ) then
      ifreq = jfreq
      exit
      endif
      enddo
      endif


c    logarithmic frequency inter/extrpolation 
      df = (Log(zfreq(ifreq+1)) - Log(zfreq(ifreq)))
      df0= (Log(freq) - Log(zfreq(ifreq)))
      FSAVE = FREQ
      endif



      THT = XTHT
      SST = XSST
      WIND = XWIND

      ITHT  = int((THT -THT0)/DTHT) + 1
      ISST  = int((SST -SST0)/DSST) + 1
      IWIND = int((WIND -WIND0)/DWIND) + 1

      if (ITHT < 1) THEN
          ITHT = 1
          THT = 0.0
      ENDIF

      if (ITHT > NTHT) THEN
          ITHT = NTHT
          THT = 88.0
      ENDIF

      if (ISST < 1) THEN
          ISST = 1
          SST = -2
      ENDIF 

      if (ISST > NSST) then
         ISST = NSST
         SST =  40
      endif

      if (IWIND < 1) THEN
          IWIND = 1
          WIND = 0
      ENDIF 

      if (IWIND > NW) then
          IWIND = NW
          WIND =  40
      ENDIF 

      THT_1  = THT0  + (ITHT - 1) *DTHT
      SST_1  = SST0  + (ISST - 1) *DSST
      WIND_1 = WIND0 + (IWIND - 1)*DWIND


      T_THT  = (THT   - THT_1)      / DTHT
      T_sst  = (sst   - sst_1)      / DSST
      T_WIND = (WIND  - WIND_1)     / DWIND

      i1 = iwind
      i2 = iwind+1
      if (i2 > NW) i2 = NW

      j1 = isst
      j2 = isst+1
      if (j2 > NSST) j2 = NSST

      k1 = itht
      k2 = itht+1
      if (k2 > ntht) k2 = ntht



      DO L = 0,1
      DO IPOL = 1,2

      FARAY(1,1) = (1-T_WIND)* EMISS_ARRAY(IFREQ+L,IPOL,k1,j1,i1) + 
     &      T_WIND * EMISS_ARRAY(IFREQ+L,IPOL,k1,j1,i2)  

      FARAY(1,2) = (1-T_WIND)* EMISS_ARRAY(IFREQ+L,IPOL,k1,j2,i1) + 
     &      T_WIND * EMISS_ARRAY(IFREQ+L,IPOL,k1,j2,i2) 

      FARAY(2,1) = (1-T_WIND)* EMISS_ARRAY(IFREQ+L,IPOL,k2,j1,i1) + 
     &      T_WIND * EMISS_ARRAY(IFREQ+L,IPOL,k2,j1,i2) 

      FARAY(2,2) = (1-T_WIND)* EMISS_ARRAY(IFREQ+L,IPOL,k2,j2,i1) + 
     &      T_WIND * EMISS_ARRAY(IFREQ+L,IPOL,k2,j2,i2) 


      GARAY(1) = (1-T_SST)*FARAY(1,1) + T_SST*FARAY(1,2)
      GARAY(2) = (1-T_SST)*FARAY(2,1) + T_SST*FARAY(2,2)

      YY(IPOL,L) = (1-T_THT)*GARAY(1) + T_THT*GARAY(2)

      ENDDO ! ipol loop
      ENDDO ! L (IFREQ , IFREQ+1) LOOP


c    FREQUENCY INTER/EXTRA POLATION
      DEMISS(1:2) = YY(1:2,0)*(1.0-DF0/DF) + YY(1:2,1)*DF0/DF
      if (DEMISS(1) < 0.) DEMISS(1) = 0.0

      END





      SUBROUTINE FDOMEGA (FREQ,THT,TRAN,WIND,  XOMEGA)
c     FDOMEGA: 
c     SCATTERTEM (non-specular reflection)
c     DEC 2002: 
c     logarithmic frequency interpolation
c     up to 240 GHZ
c     all EIA

c
c     NAME: 
c     FDOMEGA 
c     
c     USAGE:
c     CALL  (FREQ,THT,TRAN,WIND,  XOMEGA)
c
c     DESCRIPTION:
c     calculates scattertem to correct for non specular reflection for V and H pol 
c     for given EIA, FREQ, wind speed
c     see CMIS OCEAN ATBD, section 3.6 
c     logarithmic frequqency interpolation

c 
c
c
c
c    INPUT:
c    NAME       DESCRIPTION                 TYPE                   LENGTH      UNIT       RANGE
c ----------------------------------------------------------------------------------------------
c
c    FREQ        FREQUENCY                      REAL           SCALAR  GHz          [5,240]
c    THT         EIA                        REAL           SCALAR  deg      [0,88]
c    TRAN        atmospheric transmittance  REAL           SCALAR             [0,1]  
c    WIND        Ocean Surface Wind Speed   REAL           SCALAR  m/s        [0,40]
c                    10 m above ground                 
c
c    OUTPUT:
c    NAME       DESCRIPTION                           TYPE           LENGTH           RANGE
c ----------------------------------------------------------------------------------------------
c    XOMEGA     scatterterm                           REAL           VECTOR(2)            [0,1]  
c
 


      IMPLICIT NONE

      integer, parameter ::  NWIND=21,NTHT=10,NATM=20,NFREQ=9
      REAL, PARAMETER ::THT0=0.0,WIND0=0.0,TEFF=290.0,F0 = 1.5
      REAL, PARAMETER ::DTHT=10.0,DF=0.5,DWIND=1.0,DTRAN=-0.050

      REAL(4),DIMENSION(1:NFREQ,1:NWIND,1:NATM,1:NTHT,1:2),SAVE:: 
     &    SCATTERM
      ! OMEGA      ( FREQ,W,TRAN,THT,POL,T)

      INTEGER(4) :: IWIND,ITHT,ITRAN,IPOL,IFREQ

      REAL(4) :: FREQ, WIND,COSTHT,THT,PATH,TRAN,OPACTY,TRAN0
c     1PATH, TRAN, OPACTY,  TRAN0
      real(4) :: xlogf,xtht,xwind
        real    :: deg2rad

      REAL(4) :: T_THT, T_F, T_WIND,T_TRAN,F_1,THT_1, TRAN_1, WIND_1


      INTEGER :: i1,i2,j1,j2,k1,k2,l1,l2
      REAl(4), DIMENSION(2,2,2,2) :: HINT
      REAl(4), DIMENSION(2,2,2  ) :: GINT
      REAl(4), DIMENSION(2,2    ) :: FINT
      REAL(4), DIMENSION(2)       :: XOMEGA


      INTEGER, SAVE :: ISTART = 0
      SAVE


c     EXTERNAL ASCII FILE 
c     PATH + FILENAME NEED TO BE SPECIFIED BY USER     
      character (LEN=100) :: inputfile=
     &   'E:\code_for_AER\code_for_AER_2002\input_tables\scat.lis'
c     1'E:\code_for_AER\code_for_AER_2002\input_tables\scat.lis'


c     WHEN THE SUBROUTINE IS CALLED THE FIST TIME (ISTART = 0), EMISS_ARRAY 
c     are read from inputfile
c     for subsequent calls istart = 1 and the arrays from the last run are used  
      deg2rad=acos(-1.0)/180.0
      if (istart == 0) then
      open(unit=3,file=inputfile,form='formatted',status='old')
      do ifreq = 1,nfreq
      DO IPOL=1,2
      DO ITHT=1,NTHT
      do itran  = 1,natm 
      read(3,8099) SCATTERM(IFREQ,1:NWIND,ITRAN,ITHT,IPOL)
 8099 format(21(f9.4,1x))
      enddo
      enddo
      enddo
      enddo

      close(3)
      istart = 1
      endif


      XLOGF = LOG(FREQ)
      XTHT = THT
      XWIND = WIND
      costht = cos(tht*deg2rad)
      if (tran < 1.E-20) then
      Xomega = 0.0 ! opaque atmosphere
      return
      endif
      PATH=1.00035/SQRT(COSTHT*COSTHT+7.001225E-4)   
      !(1+HRATIO)/SQRT(COSTHT**2+HRATIO*(2+HRATIO)), HRATIO=.00035
      ! maximum for path is 37.80638 (if tht = 90)
      opacty = - alog(tran) /path
      tran0 = exp(-opacty)
      if (opacty < 1.0E-4) tran0 = 1.0


      IFREQ = int((XLOGF-F0)/DF)    + 1
      ITHT  = int((XTHT -THT0)/DTHT) + 1
      IWIND = int((XWIND -WIND0)/DWIND) + 1
      itran = int((tran0 - 1.0)/DTRAN) + 1


      if (ITHT < 1)      ITHT = 1
      if (ITHT > NTHT)   ITHT = NTHT

      if (ITRAN < 1)     ITRAN = 1
      if (ITRAN > NATM)  ITRAN = NATM

      if (IFREQ < 1)     IFREQ = 1
      if (IFREQ > NFREQ) IFREQ = NFREQ


      if (IWIND < 1)     IWIND = 1
      if (IWIND > NWIND) IWIND = NWIND
 
      F_1    = F0 + (IFREQ-1)*DF
      TRAN_1 = 1.0 + DTRAN*(itran-1)
      THT_1  = THT0  + (ITHT - 1) *DTHT
      WIND_1 = WIND0 + (IWIND - 1)*DWIND

      T_TRAN =  (TRAN0 - TRAN_1)  / DTRAN
      T_F    =  (XLOGF - F_1) /DF
      T_THT  = (THT   - THT_1)      / DTHT
      T_WIND = (WIND  - WIND_1)     / DWIND

      i1 = iwind
      i2 = iwind+1
      if (i2 > NWIND) i2 = NWIND

      j1 = itran
      j2 = itran+1
      if (j2 > NATM) j2 = NATM

      k1 = itht
      k2 = itht+1
      if (k2 > ntht) k2 = ntht


      l1 = ifreq
      l2 = ifreq+1
      if (l2 > nfreq) l2 = nfreq

c    wind interpolation
      HINT(1,1,1,1:2) = 
     &     (1-T_WIND)*scatterm(l1,i1,j1,k1,1:2) + 
     &     T_WIND    *scatterm(l1,i2,j1,k1,1:2)

      HINT(2,1,1,1:2) = 
     &     (1-T_WIND)*scatterm(l2,i1,j1,k1,1:2) + 
     &     T_WIND    *scatterm(l2,i2,j1,k1,1:2)

      HINT(1,1,2,1:2) = 
     &     (1-T_WIND)*scatterm(l1,i1,j1,k2,1:2) + 
     &     T_WIND    *scatterm(l1,i2,j1,k2,1:2)

      HINT(2,1,2,1:2) = 
     &     (1-T_WIND)*scatterm(l2,i1,j1,k2,1:2) + 
     &     T_WIND    *scatterm(l2,i2,j1,k2,1:2)

      HINT(1,2,1,1:2) = 
     &     (1-T_WIND)*scatterm(l1,i1,j2,k1,1:2) + 
     &     T_WIND    *scatterm(l1,i2,j2,k1,1:2)

      HINT(2,2,1,1:2) = 
     &     (1-T_WIND)*scatterm(l2,i1,j2,k1,1:2) + 
     &     T_WIND    *scatterm(l2,i2,j2,k1,1:2)

      HINT(1,2,2,1:2) = 
     &     (1-T_WIND)*scatterm(l1,i1,j2,k2,1:2) + 
     &     T_WIND    *scatterm(l1,i2,j2,k2,1:2)

      HINT(2,2,2,1:2) = 
     &     (1-T_WIND)*scatterm(l2,i1,j2,k2,1:2) + 
     &     T_WIND    *scatterm(l2,i2,j2,k2,1:2)


c    TRAN INTERPOLATION
      GINT(1,1,1:2) = 
     &     (1-T_TRAN)*HINT(1,1,1,1:2) + T_TRAN*HINT(1,2,1,1:2)
      GINT(1,2,1:2) = 
     &     (1-T_TRAN)*HINT(1,1,2,1:2) + T_TRAN*HINT(1,2,2,1:2)
      GINT(2,1,1:2) = 
     &     (1-T_TRAN)*HINT(2,1,1,1:2) + T_TRAN*HINT(2,2,1,1:2)
      GINT(2,2,1:2) = 
     &     (1-T_TRAN)*HINT(2,1,2,1:2) + T_TRAN*HINT(2,2,2,1:2)

c    THT INTERPOLATION
      FINT(1,1:2)   = (1-T_THT)*GINT(1,1,1:2) + T_THT*GINT(1,2,1:2)
      FINT(2,1:2)   = (1-T_THT)*GINT(2,1,1:2) + T_THT*GINT(2,2,1:2)

c     FREQUENCY INTERPOLATION
      XOMEGA(:)  = (1-T_F)*FINT(1,1:2) + T_F*FINT(2,1:2)

      RETURN
      END


c
      SUBROUTINE FINDTB(FREQ,THT,SST,SAL,WIND,TBUP,TBDW,TRAN,FMOD)
c     1                    FMOD)
c    RADIATIVE TRANSFER ROUTINE (CMIS OCEAN ATBD equations 3.11 + 3.62)

c    INPUT:
c    NAME       DESCRIPTION                 TYPE                   LENGTH      UNIT       RANGE
c ----------------------------------------------------------------------------------------------
c
c    FREQ        FREQUENCY                      REAL           SCALAR  GHz          [5,240]
c    THT         EIA                        REAL           SCALAR  deg      [0,88]
c    SST         Sea Surface Temperature    REAL           SCALAR  Celsius    [2,40]
c    SAL         SALINITY                   REAL           SCALAR  ppt        [0,40]           
   
c    TRAN        atmospheric transmittance  REAL           SCALAR             [0,1]  
c    WIND        Ocean Surface Wind Speed   REAL           SCALAR  m/s        [0,40]
c                    10 m above ground                 
c    TBUP        upwelling atmospheric TB   REAL           SCALAR  K 
c    TBDW      downwelling atmospheric TB   REAL           SCALAR  K 

c
c    OUTPUT:
c    NAME       DESCRIPTION                           TYPE           LENGTH          UNIT
c ----------------------------------------------------------------------------------------------
c    FMOD       RTM TB                                REAL           VECTOR(2)       Kelvin             
c               1 = v-pol,   2 = h-pol
c


      REAL(4) :: FREQ,THT,SURTEP,SST,SAL,WIND,TBUP,TBDW,TRAN,TC=2.7
        REAL(4) :: e1,e2
      REAL(4), DIMENSION(2) :: EMISS,XOMEGA,TBOMEGA,FMOD 
      SAVE

      SURTEP =  SST + 273.15
      CALL FIND_EMISS_TOT (FREQ,THT,SST,WIND,SAL,  E1,e2)
        EMISS(1)=e1
        EMISS(2)=e2
      CALL  FDOMEGA (FREQ,THT,TRAN,WIND,  XOMEGA)

      TBOMEGA = (1.0+XOMEGA)*(TBDW - (1.0-TRAN)*TC) + TC
      FMOD = TBUP + TRAN*EMISS*SURTEP + TRAN*(1.0-EMISS)*TBOMEGA


      RETURN
      end

