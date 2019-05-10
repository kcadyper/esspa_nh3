      SUBROUTINE CDMX(F1,F2,TC,NU,FLAG,K2,NAV)

C  Variables:
      implicit none
      INTEGER FLAG
      REAL F1,F2,F3,TC,K2,IUNIT,XLAM,NULL,NU
      REAL NR1,NR2,NR3,NI1,NI2,NI3,ER1,ER2,ER3,EI1,EI2,EI3,EAVR,EAVI
      REAL NAVR,NAVI,CC,ERMIX,EIMIX
      COMPLEX NAV,K
C************************************************************************
C
C  Version 2.0
C
C  This subroutine calls BRUGG or MAXGAR, EPSREFWAT, REFICE 
C  and NDX2DIE (index of refraction to dielectric constant)
C
C  It also returns the average dielectric constant, and the |K|^2  
C  parameter. It also has error checking to insure that proper values are
C  being used before they are passed into functions and subroutines.
C  DIE2NDX and NDX2DIE are subroutines used to convert back and forth
C  from dielectric constants to indices of refraction.
C 
C  Sep. 2000 : B.T. Johnson, Purdue University - Department of Earth and
C  Atmospheric Sciences.  jbenjam@eas.purdue.edu [jbenjam@aos.wisc.edu]
C
C  Revision History
C  Version 2.0  Sept. 5, 2000
C    Changed calling parameters:  Included FLAG to allow selection of 
C    mixture formula to be used.  0 = Bruggeman formula, ???0 and ???1
C    refer to Maxwell Garnett formulas, where ???0 is the average of
C    the first two components in the third component's matrix, and ???1
C    is the 1st component in the average of the 2nd and 3rd matrix.
C    e.g. 1230 = [WI]A  and 1231 = W[IA].  [WI]A is Water included in
C    Ice Matrix, which the average is then included in the Air matrix. 
C    W[IA] is Water included in the average of Ice included in Air.
C    The general trend is inclusion -> matrix from left to right.
C    
C    NOTE:  This change modifies CALLING ROUTINES because of the 
C    changes in the arguments.  
C
C  Version 1.4  August 30, 2000
C    Changed NAVR, NAVI in parameter list to complex NAV. NOTE: This
C    modifies CALLING ROUTINES because of the changes in the arguments.
C 
C  Version 1.3  August 23, 2000 
C    Added 6 additional cases to MAXGAR to account for one component
C    inclusions in a averaged matrix e.g. A[WI] = Air in Water/Ice Matrix
C    
C  Version 1.2  July 26, 2000
C    Moved subroutine into "beta" mode, corrected comments.
C  Version 1.1  February 1, 2000
C    Added Bruggeman and Maxwell Garnett subroutines.
C  Version 1.0e November 1999
C    Removed RHOAV from the passing parameters and made modificiations
C    to the notes.  Preparing for final packaging.
C  Version 1.0d August 1998
C    Included a check for the minimum and maximum allowable limits for
C    wavelength:  [0.200 - 1.0e5] microns
C    Note: this limit is only valid when water is incorporated,
C    otherwise, the limit is larger for ice: [0.045 - 8.6e6] microns.
C    The limitations are based on information from the index of
C    refraction subroutines, and do not necessarily constitute the range
C    of *physically* acceptable parameters.
C  Version 1.0c August 1998
C    Uses revised version of REFWAT to incorporate Liebe (1991) 
C    [EPSW] microwave index of refraction calculations.
C    Now valid for most wavelengths again (see Version 1.0d).
C  Version 1.0b July 1998
C    Replaced REFWAT with EPSW : Note this makes water valid only
C    for microwave wavelengths!
C    Fixed a bug in a loop: "DO I=2,3" changed to "DO I=1,3"
C  Version 1.0a July 1998
C    Original   
C 
C  The subroutine will:
C  1.  Accept inputs of F1, F2, TC, NU, FLAG. 
C  2.  Convert NU (microwave frequency) to XLAM (wavelength in 
C      in microns).  Since XLAM is in microns, IUNIT=0.  See
C      refwat.f or refice.f for more info on IUNIT and XLAM.
C  3.  Send IUNIT, XLAM, and (TC+273.15) to EPSREFWAT and REFICE to 
C      return the proper real and complex portions of the index of 
C      refraction; NR and NI.
C  4.  Convert index of refraction into the respective dielectric
C      functions using equation 9.1 on page 227 of Bohren and Huffman
C      "Absorption and Scattering of Light by Small Particles" (1983).
C      (NR --> ER,  NI --> EI) as coded in the subroutines NDX2DIE
C      and DIE2NDX.
C  5.  Send F1, F2, ERi, EIi (i = components 1 to 3 ) into BRUGG/MAXGAR,
C      and return EAVR, and EAVI.  EAVR = avg. real dielectric
C      EAVI = avg. imaginary dielectric constant.
C  6.  Convert EAVR, and EAVI into their respective indexes of refraction;
C      NAVR, and NAVI.  Combine NAVR and NAVI into one complex variable NAV,
C      yielding the complex index of refraction.
C  7.  Calculate |K|^2 (K2) using the equation K=(NAV**2-1)/(NAV**2+2) 
C      from Rogers and Yao "A Short Course in Cloud Physics".  
C      Where NAV is the average complex index of refraction of the entity. 
C  8.  Return NAV, |K|^2, to the calling procedure.
C
C  Notes:
C    > TC is degrees celcius, adding 273.15 to get Kelvins.
C    > IUNIT, XLAM  see refwat.f by Eric A. Smith or refice.f by Stephen
C      G. Warren.  XLAM is the wavelength.
C    > When compiling, be sure to include DIE2NDX, NDX2DIE, EPSREFWAT,
C      REFICE, BRUGG, and/or MAXGAR (see diemix.f)
C    > BRUGG only needs to be called once, because it mixes all three
C      components (if they exist) in one pass.  However, MAXGAR, requires
C      two calls (as is done in this routine) in order to mix the 
C      components.  The first argument supplied to MAXGAR is the inclusion 
C      volume fraction, the second argument is the matrix volume fraction.
C 
C************************************************************************
      CC=2.99792E+08
      IUNIT=0
C     NU is in GHz, so we convert to Hz and complete the division.  However,
C     XLAM needs to be in microns (1.0E-6 meters), so divide again.
      XLAM=(CC/(NU*1.0E9))*1.0E6
      if(F1.GT.0.0.AND.(XLAM.GT.1.0E6.OR.XLAM.LE.0.200))
     &write(*,*)'WAVELENGTH OUT OF RANGE:Water [0.200 um-1.0e+5 um]'
      if(F2.GT.0.0.AND.(XLAM.GT.8.6E6.OR.XLAM.LE.0.045))
     &write(*,*)'WAVELENGTH OUT OF RANGE: ice  [0.045 um-8.6e+6 um]'

      call EPSREFWAT(IUNIT,XLAM,(TC+273.15),NR1,NI1,NULL,NULL)
      call REFICE(IUNIT,XLAM,(TC+273.15),NR2,NI2,NULL,NULL)

C Dielectric constant of air... assumed to be 1.0 since the surrounding
C medium is also air, and is assumed to be 1.0 in Mie calculations.


      ER3 = 1.0
      EI3 = 0.0
      F3  = (1-F1-F2)
      call NDX2DIE(NR1,NI1,ER1,EI1)
      call NDX2DIE(NR2,NI2,ER2,EI2)
C      call NDX2DIE(NR3,NI3,ER3,EI3)
      call DIE2NDX(ER3,EI3,NR3,NI3)

C Use Bruggeman Formula
      if (flag .eq. 0) then
         call BRUGG(F1,F2,ER1,EI1,ER2,EI2,ER3,EI3,EAVR,EAVI)
      else if (flag .eq. 1230) then
C     [WI]A
         if (F1 .gt. 0.0 .or. F2. gt. 0.0) then             
            call MAXGAR(F1,F2,ER1,EI1,ER2,EI2,EAVR,EAVI)
         end if
         if ((1-F3) .lt. 1.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR((1-F3),F3,ERMIX,EIMIX,ER3,EI3,EAVR,EAVI)
         end if
      else if (flag .eq. 2130) then
C     [IW]A
         if (F1 .gt. 0.0 .or. F2. gt. 0.0) then      
            call MAXGAR(F2,F1,ER2,EI2,ER1,EI1,EAVR,EAVI)
         end if
         if ((1-F3) .lt. 1.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR((1-F3),F3,ERMIX,EIMIX,ER3,EI3,EAVR,EAVI)
         end if
      else if (flag .eq. 3210) then
C     [AI]W
         if (F3 .gt. 0.0 .or. F2. gt. 0.0) then
            call MAXGAR(F3,F2,ER3,EI3,ER2,EI2,EAVR,EAVI)
         end if
         if (F1 .gt. 0.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR((F2+F3),F1,ERMIX,EIMIX,ER1,EI1,EAVR,EAVI)
         end if
      else if (flag .eq. 2310) then         
C     [IA]W
         if (F3 .gt. 0.0 .or. F2. gt. 0.0) then       
            call MAXGAR(F2,F3,ER2,EI2,ER3,EI3,EAVR,EAVI)
         end if
         if (F1 .gt. 0.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR((F2+F3),F1,ERMIX,EIMIX,ER1,EI1,EAVR,EAVI)
         end if
      else if (flag .eq. 1320) then
C     [WA]I
         if (F1 .gt. 0.0 .or. F3. gt. 0.0) then
            call MAXGAR(F1,F3,ER1,EI1,ER3,EI3,EAVR,EAVI)
         end if
         if (F2 .gt. 0.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR((F1+F3),F2,ERMIX,EIMIX,ER2,EI2,EAVR,EAVI)
         end if
      else if (flag .eq. 3120) then
C     [AW]I
         if (F3 .gt. 0.0 .or. F1. gt. 0.0) then
            call MAXGAR(F3,F1,ER3,EI3,ER1,EI1,EAVR,EAVI)
         end if
         if (F2 .gt. 0.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR((F1+F3),F2,ERMIX,EIMIX,ER2,EI2,EAVR,EAVI)
         end if
      else if (flag .eq. 1231) then
C     W[IA]
         if (F2 .gt. 0.0 .or. F3. gt. 0.0) then
            call MAXGAR(F2,F3,ER2,EI2,ER3,EI3,EAVR,EAVI)
         end if
         if (F1 .gt. 0.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR(F1,(F2+F3),ER1,EI1,ERMIX,EIMIX,EAVR,EAVI)
         end if
      else if (flag .eq. 2131) then
C     I[WA]
         if (F3 .gt. 0.0 .or. F1. gt. 0.0) then
            call MAXGAR(F1,F3,ER1,EI1,ER3,EI3,EAVR,EAVI) 
         end if
         if (F2 .gt. 0.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR(F2,(F1+F3),ER2,EI2,ERMIX,EIMIX,EAVR,EAVI)
         end if
      else if (flag .eq. 3211) then
C     A[IW]
         if (F2 .gt. 0.0 .or. F1. gt. 0.0) then
            call MAXGAR(F2,F1,ER2,EI2,ER1,EI1,EAVR,EAVI)
         end if
         if ((1.0-F3) .lt. 1.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR(F3,(1.0-F3),ER3,EI3,ERMIX,EIMIX,EAVR,EAVI)
         end if
      else if (flag .eq. 2311) then
C     I[AW]
         if (F3 .gt. 0.0 .or. F1. gt. 0.0) then
            call MAXGAR(F3,F1,ER3,EI3,ER1,EI1,EAVR,EAVI)
         end if
         if (F2 .gt. 0.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR(F2,(F1+F3),ER2,EI2,ERMIX,EIMIX,EAVR,EAVI)
         end if
      else if (flag .eq. 1321) then
C     W[AI]
         if (F3 .gt. 0.0 .or. F2. gt. 0.0) then
            call MAXGAR(F3,F2,ER3,EI3,ER2,EI2,EAVR,EAVI)
         end if
         if (F1 .gt. 0.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR(F1,(F2+F3),ER1,EI1,ERMIX,EIMIX,EAVR,EAVI)
         end if
      else if (flag .eq. 3121) then         
C     A[WI]
         if (F2 .gt. 0.0 .or. F1. gt. 0.0) then
            call MAXGAR(F1,F2,ER1,EI1,ER2,EI2,EAVR,EAVI)
         end if
         if ((1.0-F3) .lt. 1.0) then
            ERMIX = EAVR
            EIMIX = EAVI
            call MAXGAR(F3,(1.0-F3),ER3,EI3,ERMIX,EIMIX,EAVR,EAVI)
         end if
      end if	

      call DIE2NDX(EAVR,EAVI,NAVR,NAVI)
      NAV=CMPLX(NAVR,NAVI)

      K=(NAV*NAV-1)/(NAV*NAV+2)
      K2=(cabs(K))**2
      return 
      end







