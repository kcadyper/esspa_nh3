      SUBROUTINE BRUGG(F1,F2,ER1,EI1,ER2,EI2,ER3,EI3,EAVR,EAVI)
C  Variables:
      IMPLICIT none
      REAL F1,F2,F
      REAL ER1,ER2,ER3,EI1,EI2,EI3,EAVR,EAVI
      COMPLEX E1,E2,E3,EAV,B,EAV0
C*************************************************************************
C  Subroutine calculates the average dielectric function for up to 3 
C  number of constituents.  It uses a modified version of the Bruggeman 
C  dielectric mixing function. (modified to use 3 materials instead of 2)
C  Called by CDMX subroutine.
C
C  July 2000 : Benjamin T. Johnson, Atmospheric and Oceanic Sciences Dept.
C              University of Wisconsin - Madison 
C              jbenjam@aos.wisc.edu
C 
C  The subroutine will:
C  1. Accept the parameters F1, F2, ER1,EI1,ER2,EI2,ER3,EI3.   
C
C  2. Calculate the real and complex portions of the average or mixed
C     dielectric constant, and return it back to the calling program.
C
C  Variable and parameter descriptions:  
C     F1 = volume fraction of component 1
C     F2 = volume fraction of component 2
C     1-F1-F2 = volume fraction of component 3 
C     F  = volume fraction of inclusion.
C     
C     ERx, EIx are the real and complex parts of the dielectric  
C     constant for each of the three materials (x=1,2,3)
C 
C     EAVR, EAVI are average or mixed real and imaginary parts of the 
C     dielectric constant.
C*************************************************************************
      F=F1/(F1+F2)
      E1=cmplx(ER1,EI1)
      E2=cmplx(ER2,EI2)
      E3=cmplx(ER3,EI3)
      B=((1.0-3.0*(1.0-F))*E2+(1.0-3.0*F)*E1)/2.0
      EAV0=-0.5*B+0.5*(B*B+2.0*E2*E1)**0.5
      F=(F1+F2)/1.0
      B=((1.0-3.0*(1.0-F))*E3+(1.0-3.0*F)*EAV0)/2.0
      EAV=-0.5*B+0.5*(B*B+2.0*E3*EAV0)**0.5
C Converts back to real and imaginary components
c      EAVR=sngl(EAV)
      EAVR=real(EAV)
      EAVI=aimag(EAV)
      return
      end

C ========================================================================
C ========================================================================

      SUBROUTINE MAXGAR(FI,FM,ERI,EII,ERM,EIM,EAVR,EAVI)
C  Variables:
      IMPLICIT none
      REAL FI,FM,F
      REAL ERI,ERM,EII,EIM,EAVR,EAVI
      COMPLEX EI,EM,EAV,B
C*************************************************************************
C  Subroutine calculates the average dielectric function for TWO 
C  number of constituents.  It is based on the Maxwell Garnett 
C  dielectric mixing function. >>>ONLY MIXES TWO MATERIALS<<<<
C
C  Note: When mixing materials, considerations must be made as to which
C        component will be the inclusion or the matrix.  Significant 
C        differences appear with respect to the average dielectric
C        if the two materials are simply interchanged.  By this same
C        token, care must be taken when generalizing to more than two
C        components as to which materals will be matrix and inclusion.
C        See Bohren and Huffman, 1983: Absorption and Scattering of Light
C        by small particles.  The Bruggeman routine (BRUGG) does not have
C        the distinction of inclusion and matrix. 
C
C  July 2000 : Benjamin T. Johnson [jbenjam@aos.wisc.edu] 
C  
C     INPUTS: FI, FM, ERI,EII,ERM,EIM
C     OUTPUT: EAVR, EAVI
C
C  The subroutine will:
C  1. Accept the parameters FI,FM, ERI,ERM,EII,EIM  
C
C  2. Calculate the real and complex portions of the average or mixed
C     dielectric constant, and return it back to the calling program.
C
C  Variable and parameter descriptions:  
C     FI is the fractional part of component that is the inclusion
C     with respect to component FM. FM is the fraction of components
C     that make up the "matrix" of which FI is included.
C     
C     ERI, EII are the real and complex parts of the dielectric  
C     constant for the included material, and ERM, EIM are for the
C     matrix material.
C 
C     EAVR, EAVI are average or mixed real and imaginary parts of the 
C     dielectric constant.
C 
C*************************************************************************
C F is the volume fraction of the inclusions, thus, FI is the inclusion.
      F    = FI/(FI+FM) 
      EI   = cmplx(ERI,EII)
      EM   = cmplx(ERM,EIM)
      B    = ((EI-EM)/(EI+2.0*EM))
      EAV  = EM*(1.0+(3.0*F*B)/(1.0-F*B))
c      EAVR = sngl(EAV)
      EAVR = real(EAV)
      EAVI = aimag(EAV)
      return
      end

