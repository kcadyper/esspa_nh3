      subroutine rosenk(t,p,rho,f,attn,absh2o_l,absh2o_c,absorpN2,
     &     absO2)
      REAL  t,p,rho,f,attn,absh2o_l,absh2o_c,absorpN2,absO2
      !------------------------------------------------
      !Rosenkranz subroutines reorganized to accept 
      !Basic water vapor density (not mixing ratio).
      !Minor changes were also made in order to use
      !more precise constants (217 and 0.5034) to be
      !compatible with the subroutines that accept
      !mass amount as input.
      !S.A.B. 2003 AER Inc.
      !------------------------------------------------
      !Inputs:
      !-------
      !t   : temperature in K
      !p   : total pressure in mb
      !rho : water vapor density in g/m3
      !f   : frequency in GHz
      !
      !Outputs:
      !--------
      !attn     : total absorption [Nepers/km]
      !absh2o_l : Water Vapor absorption due to lines [Nepers/km]
      !absh2o_c : Water Vapor absorption due to cont. [Nepers/km]
      !absorpN2 : N2-induced absorption [Nepers/km]
      !absO2    : O2-induced absorption [Nepers/km]
      !------------------------------------------------
      absh2o_l  = ABH2O_l(t,p,rho,f)
      absh2o_c  = ABH2O_c(t,p,rho,f)
      absorpN2  = ABSN2(t,p,f)
      absO2     = O2ABS(t,p,rho,f)
      attn      = absh2o_l+absh2o_c+absorpN2+absO2
      return
      end


      subroutine rosenk_u(t,p,uh2o,udry,o2ratio,f,attn,absh2o_l,
     &     absh2o_c,absorpN2,absO2)
      implicit none
      REAL t,p,uh2o,udry,f,attn,absh2o_l,absh2o_c,absorpN2,absO2
      REAL ABH2O_l_u,ABH2O_c_u,ABSN2_u,O2ABS_u,o2ratio
      !------------------------------------------------
      !Variant of Rosenkranz subroutines reorganized to 
      !accept water vapor and dry air mass amounts.
      !See notebook for formalism used for derivation
      !S.A.B. 2003 AER Inc.
      !This was developed to be compatible with the 
      !new RTM scheme
      !------------------------------------------------
      !Inputs:
      !-------
      !t    : temperature in K
      !p    : total pressure in mb
      !uh2o : water vapor mass amount in g/cm2
      !udry : dry air mass amount in g/cm2
      !f    : frequency in GHz
      !O2Ratio : O2 ratio wrt dry air (replaces the default 0.2085) 
      !
      !Outputs:
      !--------
      !attn     : total optical depth [neper]
      !absh2o_l : Water Vapor optical depth due to lines [neper]
      !absh2o_c : Water Vapor optical depth due to cont. [neper]
      !absorpN2 : N2-induced optical depth [neper]
      !absO2    : O2-induced optical depth [neper]
      !------------------------------------------------
      absh2o_l  = ABH2O_l_u(t,p,uh2o,udry,f)
      absh2o_c  = ABH2O_c_u(t,p,uh2o,udry,f)
      absorpN2  = ABSN2_u(t,p,uh2o,udry,f)
      absO2     = O2ABS_u(t,p,uh2o,udry,o2ratio,f)
      attn      = absh2o_l+absh2o_c+absorpN2+absO2
      return
      end





      FUNCTION ABH2O_l(T,P,RHO,F)
      IMPLICIT NONE
      REAL T,P,RHO,F,ABH2O_l	
C      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
C      T       KELVIN    I   TEMPERATURE
C      P       MILLIBAR  I   PRESSURE              .1 TO 1000
C      RHO     G/M**3    I   WATER VAPOR DENSITY
C      F       GHZ       I   FREQUENCY             0 TO 800
C      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
C   LOCAL VARIABLES:
      INTEGER NLINES,I,J
      PARAMETER (NLINES=15)
      REAL DF(2),S1(NLINES),B2(NLINES),W3(NLINES),FL(NLINES),X(NLINES),
     & WS(NLINES),XS(NLINES)
      REAL PVAP,PDA,DEN,TI,TI2,SUM,WIDTH,WSQ,S,BASE,RES,CON
C     LINE FREQUENCIES:
      DATA FL/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,
     & 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,
     & 620.7008, 752.0332, 916.1712/
C     LINE INTENSITIES AT 300K:
      DATA S1/ .1310E-13, .2273E-11, .8036E-13, .2694E-11, .2438E-10,
     & .2179E-11, .4624E-12, .2562E-10, .8369E-12, .3263E-11, .6659E-12,
     & .1531E-08, .1707E-10, .1011E-08, .4227E-10/
C     T COEFF. OF INTENSITIES:
      DATA B2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,
     & 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
C     AIR-BROADENED WIDTH PARAMETERS AT 300K:
      DATA W3/.00281, .00281, .0023, .00278, .00287, .0021, .00186,
     & .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
C     T-EXPONENT OF AIR-BROADENING:
      DATA X/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69,
     & .71, .68, .70/
C     SELF-BROADENED WIDTH PARAMETERS AT 300K:
      DATA WS/.01349, .01491, .0108, .0135, .01541, .0090, .00788,
     & .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
C     T-EXPONENT OF SELF-BROADENING:
      DATA XS/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72,
     & 1.0, .68, .84, .78/
C
      IF(RHO.LE.0.) THEN
        ABH2O_l = 0.
        RETURN
      ENDIF
      !----test-----
      !PVAP = RHO*T/217.
      PVAP = RHO*T/216.68243 !217. was an approximation
      !-------------
      PDA = P -PVAP
      DEN = 3.335E16*RHO
      TI = 300./T
      TI2 = TI**2.5
C     ADD RESONANCES
      SUM = 0.
      DO 30 I=1,NLINES
      WIDTH = W3(I)*PDA*TI**X(I) + WS(I)*PVAP*TI**XS(I)
      WSQ = WIDTH*WIDTH
      S = S1(I)*TI2*EXP(B2(I)*(1.-TI))
      DF(1) = F - FL(I)
      DF(2) = F + FL(I)
C  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
      BASE = WIDTH/(562500. + WSQ)
C  DO FOR POSITIVE AND NEGATIVE RESONANCES
      RES = 0.
      DO 20 J=1,2
      IF(ABS(DF(J)).LT.750.) RES = RES + WIDTH/(DF(J)**2+WSQ) - BASE
20    CONTINUE
30    SUM = SUM + S*RES*(F/FL(I))**2
      ABH2O_l = .3183E-4*DEN*SUM 
      RETURN
      END


      FUNCTION ABH2O_l_u(T,P,UH2O,UDRY,F)
      IMPLICIT NONE
      REAL T,P,UH2O,UDRY,F,ABH2O_l_u
      REAL H2OMIXRATIO,WVMOLMASS,DRYMOLMASS
C      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
C      T       KELVIN    I   TEMPERATURE
C      P       MILLIBAR  I   PRESSURE              .1 TO 1000
C      UH2O    G/CM**2   I   WATER VAPOR MASS AMOUNT
C      UDRY    G/CM**2   I   DRY AIR MASS AMOUNT
C      F       GHZ       I   FREQUENCY             0 TO 800
C      ABH2O   NEPERS    O   OPTICAL DEPTH
C   LOCAL VARIABLES:
      INTEGER NLINES,I,J
      PARAMETER (NLINES=15)
      REAL DF(2),S1(NLINES),B2(NLINES),W3(NLINES),FL(NLINES),X(NLINES),
     & WS(NLINES),XS(NLINES)
      REAL PVAP,PDA,DEN,TI,TI2,SUM,WIDTH,WSQ,S,BASE,RES,CON
C     LINE FREQUENCIES:
      DATA FL/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,
     & 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,
     & 620.7008, 752.0332, 916.1712/
C     LINE INTENSITIES AT 300K:
      DATA S1/ .1310E-13, .2273E-11, .8036E-13, .2694E-11, .2438E-10,
     & .2179E-11, .4624E-12, .2562E-10, .8369E-12, .3263E-11, .6659E-12,
     & .1531E-08, .1707E-10, .1011E-08, .4227E-10/
C     T COEFF. OF INTENSITIES:
      DATA B2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,
     & 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
C     AIR-BROADENED WIDTH PARAMETERS AT 300K:
      DATA W3/.00281, .00281, .0023, .00278, .00287, .0021, .00186,
     & .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
C     T-EXPONENT OF AIR-BROADENING:
      DATA X/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69,
     & .71, .68, .70/
C     SELF-BROADENED WIDTH PARAMETERS AT 300K:
      DATA WS/.01349, .01491, .0108, .0135, .01541, .0090, .00788,
     & .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
C     T-EXPONENT OF SELF-BROADENING:
      DATA XS/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72,
     & 1.0, .68, .84, .78/
      DATA WVMOLMASS /18.016 /, DRYMOLMASS/28.97/   

      IF(UH2O.LE.0.) THEN
        ABH2O_l_u = 0.
        RETURN
      ENDIF
      H2OMIXRATIO=UH2O/UDRY
      PVAP = P*(H2OMIXRATIO/(H2OMIXRATIO+(WVMOLMASS/DRYMOLMASS)))
      PDA = P -PVAP
      DEN = 3.335E16*UH2O
      TI = 300./T
      TI2 = TI**2.5
      !---ADD RESONANCES
      SUM = 0.
      DO I=1,NLINES
         WIDTH = W3(I)*PDA*TI**X(I) + WS(I)*PVAP*TI**XS(I)
         WSQ = WIDTH*WIDTH
         S = S1(I)*TI2*EXP(B2(I)*(1.-TI))
         DF(1) = F - FL(I)
         DF(2) = F + FL(I)
         !---USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
         BASE = WIDTH/(562500. + WSQ)
         RES = 0.               
         DO J=1,2
            IF(ABS(DF(J)).LT.750.) RES=RES+WIDTH/(DF(J)**2+WSQ)-BASE
         ENDDO
         SUM = SUM + S*RES*(F/FL(I))**2
      ENDDO
      ABH2O_l_u = .3183E-3*DEN*SUM 
      RETURN
      END




      FUNCTION ABH2O_c(T,P,RHO,F)
      IMPLICIT NONE
      REAL T,P,RHO,F,ABH2O_c,PVAP,PDA,TI,CON
      INTEGER I,J
      IF(RHO.LE.0.) THEN
        ABH2O_c = 0.
        RETURN
      ENDIF
      !-----test----
      !PVAP = RHO*T/217.
      PVAP = RHO*T/216.68243 !217 was an approximation: 
      !-------------
      PDA = P -PVAP
      TI = 300./T
      CON = (5.43E-10*PDA*TI**3 + 1.8E-8*PVAP*TI**7.5) 
      ABH2O_c=CON*PVAP*F*F
      RETURN
      END


      FUNCTION ABH2O_c_u(T,P,UH2O,UDRY,F)
      IMPLICIT NONE
      REAL T,P,UH2O,UDRY,F,ABH2O_c_u,PVAP,PDA,TI,CON
      REAL H2OMIXRATIO,WVMOLMASS,DRYMOLMASS
      INTEGER I,J
      DATA WVMOLMASS /18.016 /, DRYMOLMASS/28.97/   

      IF(UH2O.LE.0.) THEN
        ABH2O_c_u = 0.
        RETURN
      ENDIF
      H2OMIXRATIO=UH2O/UDRY
      PVAP = P*(H2OMIXRATIO/(H2OMIXRATIO+(WVMOLMASS/DRYMOLMASS)))
      PDA = P -PVAP
      TI = 300./T
      CON=(5.43E-10*PDA*TI**2 + 1.8E-8*PVAP*TI**6.5)
      ABH2O_c_u=CON*13.845146*UH2O*F*F !1.3845..=300/217
      RETURN
      END


      FUNCTION ABSN2(T,P,F)
C     ABSN2 = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR
C             (NEPER/KM)
C     T = TEMPERATURE (K)
C     P = PRESSURE (MB)
C     F = FREQUENCY (GHZ)
C
      real f,absn2
      TH = 300./T
      X=6.4E-14*P*F*F*TH**2.55
      ABSN2 = X*P*TH
      RETURN
      END

      FUNCTION ABSN2_u(T,P,UH2O,UDRY,F)
      Implicit none
C     ABSN2 = OPTICAL DEPTH DUE TO NITROGEN IN AIR(NEPER)
C     T = TEMPERATURE (K)
C     P = PRESSURE (MB)
C     F = FREQUENCY (GHZ)
C     UH2O    G/CM**2   WATER VAPOR MASS AMOUNT
C     UDRY    G/CM**2   DRY AIR MASS AMOUNT
C
      real t,p,f,absn2_u,UH2O,UDRY,th,WVMOLMASS,DRYMOLMASS,xKB,xNa,x
      DATA WVMOLMASS /18.016 /, DRYMOLMASS/28.97/
      DATA xKB/1.3806503E-16/, xNa/6.02214199E+23/

      TH = 300./T
      X       = 6.4E-14*P*F*F*TH**2.55
      ABSN2_u =(X*3.E-6)*(xkB*xNa)*((UDRY/DRYMOLMASS)+(UH2O/WVMOLMASS))
      RETURN
      END

      BLOCK DATA
      COMMON /O2COM/ X,WB300,W300(40),F(40),Y300(40),S300(40),
     & V(40),BE(40)
C      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      DATA F/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,
     2  59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
     3  56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,
     4  55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,
     5  53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,
     6  52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7631,
     7  487.2494, 715.3932, 773.8397, 834.1453/
        DATA S300/.2936E-14,.8079E-15, .2480E-14,.2228E-14,
     &  .3351E-14,.3292E-14, .3721E-14,.3891E-14,
     &  .3640E-14,.4005E-14, .3227E-14,.3715E-14,
     &  .2627E-14,.3156E-14, .1982E-14,.2477E-14,
     &  .1391E-14,.1808E-14, .9124E-15,.1230E-14,
     &  .5603E-15,.7842E-15, .3228E-15,.4689E-15,
     &  .1748E-15,.2632E-15, .8898E-16,.1389E-15,
     &  .4264E-16,.6899E-16, .1924E-16,.3229E-16,
     &  .8191E-17,.1423E-16, .6460E-15, .7047E-14, .3011E-14,
     &  .1826E-14, .1152E-13, .3971E-14/
      DATA BE/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626,
     & 2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625,
     & 2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844,
     & 2*7.744, .048, .044, .049, .145, .141, .145/
C      WIDTHS IN MHZ/MB
      DATA WB300/.56/, X/.8/
      DATA W300/1.63, 1.646, 1.468, 1.449, 1.382, 1.360,
     & 1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,
     & 1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05,
     & 2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.92, 3*1.81/
      DATA Y300/  -0.0233,  0.2408, -0.3486,  0.5227,
     & -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,
     &  0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,
     &  0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,
     &  0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,
     &  0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
      DATA V/  0.0079, -0.0978,  0.0844, -0.1273,
     &  0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,
     &  0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,
     &  0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,
     &  0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,
     &  0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./
      END BLOCK DATA
C
      FUNCTION O2ABS(TEMP,PRES,VAPDEN,FREQ)

C
C     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
C              IN NEPERS/KM
C
C      5/1/95  P. Rosenkranz
C
C     ARGUMENTS:
      REAL TEMP,PRES,VAPDEN
C
C     NAME    UNITS    DESCRIPTION        VALID RANGE
C
C     TEMP    KELVIN   TEMPERATURE        (UNCERTAIN)
C     PRES   MILLIBARS PRESSURE           (3 TO 1000)
C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
C     FREQ    GHZ      FREQUENCY          (0 TO 900)
C
C     REFERENCE FOR EQUATIONS AND COEFFICIENTS:
C     P.W. ROSENKRANZ, CHAP. 2 AND APPENDIX, IN ATMOSPHERIC REMOTE SENSING
C      BY MICROWAVE RADIOMETRY (M.A. JANSSEN, ED. 1993)
C     AND H.J. LIEBE ET AL, JQSRT V.48, PP.629-643 (1992)
C     (EXCEPT: SUBMILLIMETER LINE INTENSITIES FROM HITRAN92)
C
      TH = 300./TEMP
      TH1 = TH-1.
      B = TH**X
      !----test----
      !PRESWV = VAPDEN*TEMP/217.
      PRESWV = VAPDEN*TEMP/216.68243 !217 was an approximation
      !------------
      PRESDA = PRES -PRESWV
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
      DFNR = WB300*DEN
      SUM = 1.6E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
      DO 32 K=1,40
      DF = W300(K)*DEN
      Y = .001*PRES*B*(Y300(K)+V(K)*TH1)
      STR = S300(K)*EXP(-BE(K)*TH1)
      SF1 = (DF + (FREQ-F(K))*Y)/((FREQ-F(K))**2 + DF*DF)
      SF2 = (DF - (FREQ+F(K))*Y)/((FREQ+F(K))**2 + DF*DF)
32    SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))**2
      !----test----
      !O2ABS = .5034E12*SUM*PRESDA*TH**3/3.14159
      O2ABS = .503386E12*SUM*PRESDA*TH**3/3.14159 !.5034 replaced by 0.2085/(300.*xkB)
      !-----------
      RETURN
      END


      FUNCTION O2ABS_u(TEMP,PRES,UH2O,UDRY,O2Ratio,FREQ)
C     RETURNS OPTICAL DEPTH DUE TO OXYGEN IN AIR, IN NEPERS
      REAL TEMP,PRES,UH2O,UDRY,FREQ,O2ABS_u
      REAL WVMOLMASS, DRYMOLMASS, xNa,O2Ratio
C     NAME    UNITS    DESCRIPTION        VALID RANGE
C     TEMP    KELVIN   TEMPERATURE        (UNCERTAIN)
C     PRES   MILLIBARS PRESSURE           (3 TO 1000)
C     UH2O    G/CM**2   WATER VAPOR MASS AMOUNT
C     UDRY    G/CM**2   DRY AIR MASS AMOUNT
C     O2Ratio Mol/Mol   O2 ratio wrt dry air (replaces the default 0.2085) 
C     FREQ    GHZ      FREQUENCY          (0 TO 900)

      DATA WVMOLMASS /18.016 /, DRYMOLMASS/28.97/, xNa/6.02214199E+23/   
C
      TH = 300./TEMP
      TH1 = TH-1.
      B = TH**X
      H2OMIXRATIO=UH2O/UDRY
      PRESWV=PRES*(H2OMIXRATIO/(H2OMIXRATIO+(WVMOLMASS/DRYMOLMASS)))
      PRESDA = PRES -PRESWV
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
      DFNR = WB300*DEN
      SUM = 1.6E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
      DO K=1,40
         DF = W300(K)*DEN
         Y = .001*PRES*B*(Y300(K)+V(K)*TH1)
         STR = S300(K)*EXP(-BE(K)*TH1)
         SF1 = (DF + (FREQ-F(K))*Y)/((FREQ-F(K))**2 + DF*DF)
         SF2 = (DF - (FREQ+F(K))*Y)/((FREQ+F(K))**2 + DF*DF)
         SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))**2
      ENDDO
      !O2Ratio = .2085
      O2ABS_u = O2Ratio*1.E-9*SUM*(UDRY/DRYMOLMASS)*xNa*TH**2/3.14159
      RETURN
      END


