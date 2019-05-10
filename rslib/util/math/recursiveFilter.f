      SUBROUTINE RFAN1Z(AM,NAI,NSI,NX,NY,NOB  
     &     ,DELX,DELY,RM_IN,T0,T_IN,RF,WBG,XOB,YOB,OB,WOB,
     &      AA,AI,W,B1,B2,AO,QOB)     
c***********************************************************************
c* Function Name: rfan1z
c* Purpose: Run the 2-dimensional recursive filter analysis
c* Description:
c*     Formulas and methods used -
c*        All the formulas are as described in the reference, with the 
c*        following exceptions,
c*        1.) Eq(8) is modified to read:
c*               An+1 = An + (C*QW(O-An) + Wb(Ab-An)) / (C*QW + Wb) ,
c*            where Ab is the background, An the analysis after the nth
c*            increment, Wb the weight assigned to the background.
c*        2.) The radius of influence R and tolerance T are computed 
c*            for the nth increment using prescribed scaling factors.
c*            (default values are from Table C1 in HayP95)
c*        3.) In Eq(B5) in the paper, Rf*Rn should be replaced by
c*            delta*Rf*Rn.
c*
c* Inputs:
c*  Var_name          Description
c*  --------          -----------
c* AM(NX,NY)          First guess and background field.
c* NAI                Number of successive increments.
c* NSI                Number of passes of the 2-D filter (for each increment).
c* NX, NY             Number of grid points in x, y
c* NOB                Number of observations
c* DELX(NY),DELY(NY)  Grid spacing in x, y directions (assumed equal): 
c*                    dely(j)=y(j+1)-y(j), dely(ny)=dely(ny-1)
c* XOB(NOB),YOB(NOB)  Normalized x, y coordinates of
c*                    observations, x=xi + (xob-i)*DELX(j),
c*                    y=yj+(yob-j)*DELY(j), where (i,j) = int(xob,yob)+1 and
c*                    (xi,yj) are the x, y coordinates of the analysis grid
c*                    point (i,j).  XOB and YOB are only used for the bilinear
c*                    interpolation, using the 4 points (i,j), (i+1,j),
c*                    (i,j+1), (i+1,j+1).
c*                    Observations with xob<0 will be skipped.
c* OB(NOB),WOB(NOB)   Observation values and relative weights.
c* RM_IN              Scaling factors for maximum radius of influence 
c*                    (in units of DELX, DELY) for each increment.
c* T0                 Initial tolerance for residuals (in same units as OB).
c* T_IN(NAI)          Scaling factors for tolerance (relative to T0).
c* RF                 Factor for sharpening or broadening radius of influence.
c*                    Suggested value- 2 for noisy data, 0.5 for clean data.
c*                    Note- Rmax=R0/RF, Rlow=RL/RF
c* Outputs:
c* AA(NX,NY)          Final analysis field
c* AI(NX,XY)          (Work space) On output, it contains the filtered 
c*                    product of weight and residual used in the last 
c*                    increment
c* W(NX,XY)           (Work space) On output, it contains the filtered
c*                    weights used in the last increment
c* B1(NX,NY)          (Work space) On output, it contains the values of 
c*                    b=(1-alfa) for the x-direction used in the last 
c*                    increment
c* B2(NX,NY)          (Work space) On output, it contains the values of
c*                    b=(1-alfa) for the y-direction used in the last 
c*                    increment
c* AO(NOB)            (Work space) On output, it contains the f.g. values
c*                    at the observation locations used in the last increment
c* QOB(NOB)           (Work space) On output, it contains QW, the product
c*                    of obs quality and relative weight, at the observation
c*                    locations, used in the last increment
c***********************************************************************
C  ...All 2-D arrays are treated as 1-D arrays locally, and array 
C  ...indices computed accordingly
      real  AA(NX*NY),AI(NX*NY),W(NX*NY),B1(NX*NY),B2(NX*NY)
     &     ,AM(NX*NY),AO(NOB),QOB(NOB),XOB(NOB),YOB(NOB),OB(NOB) 
     &     ,WOB(NOB),WMA(NAI),RM_IN(NAI),T_IN(NAI)
      NXY=NX*NY  
      R0=RM_IN(1)*delx
C  ...Initialize analysis to background, weights to dummy value:
      DO 203 LXY=1,NXY  
         AA(LXY)=AM(LXY) 
         W(LXY)=1000.  
  203 CONTINUE
C  ...Loop over analysis increments:
      DO 300 IAI=1,NAI  
         RM=RM_IN(IAI)*delx
         WMA(iai)=DELX*DELY/(RM*RM) ! Weight given to background for each increment
         T=T0*T_IN(IAI)
         do 206 j=1,ny
C        ...Precompute constants for smoothing parameter b=1-alfa:
            DELXY=DELX*DELY
            R2L1=SQRT(2.*NSI)*DELX
            R2L2=SQRT(2.*NSI)*DELY
            RLOF1=.5*(R2L1/RF)**2/DELXY  
            RLOF2=.5*(R2L2/RF)**2/DELXY  
C        ...Eq(B4):
            E1=.5*(R2L1/R0)**2  
            E2=.5*(R2L2/R0)**2  
C        ...Eq(B3):
            B10=SQRT(E1*(E1+2.))-E1  
            B20=SQRT(E2*(E2+2.))-E2  
            WMIN=DELXY*(RF/R0)**2  
C        ...Compute smoothing parameter b=1-alfa:
C        ...Eq(B5):
            WMAX=DELXY*(RF/RM)**2  
            E1=.5*(R2L1/RM)**2  
            E2=.5*(R2L2/RM)**2  
C        ...Eq(B3):
            B1S=SQRT(E1*(E1+2.))-E1  
            B2S=SQRT(E2*(E2+2.))-E2  
C        ...Loop over grid points:
            DO 205 i=1,nx
               LXY=(j-1)*nx + i
               IF(W(LXY).GT.WMAX)THEN 
C              ...Use R=RM for very data dense areas:
                  B1(LXY)=B1S  
                  B2(LXY)=B2S  
               ELSEIF(W(LXY).LT.WMIN)THEN
C              ...Use R=R0 for very data sparse areas:
                  B1(LXY)=B10  
                  B2(LXY)=B20  
               ELSE 
C              ...Use R=sqrt(DELX*DELY)*RF*RN for the rest:
C              ...Eq(B4):
                  WSCALE=W(LXY) 
                  E1=RLOF1*WSCALE  
                  E2=RLOF2*WSCALE  
C              ...Eq(B3):
                  B1(LXY)=SQRT(E1*(E1+2.))-E1  
                  B2(LXY)=SQRT(E2*(E2+2.))-E2  
               ENDIF  
C           ...Initialize increments and weights to zero:
               W(LXY)=0.  
               AI(LXY)=0.  
C        ...Endloop over grid points:
  205       CONTINUE  
  206    CONTINUE
C     ...Compute f.g. values at obs locations:
         CALL GTOO(XOB,YOB,NOB,NX,NY,AA,AO)
         DO 251 JOB=1,NOB  
            if(xob(job).le.0) goto 251
            ODIF=OB(JOB)-AO(JOB) 
C        ...Eq(7):
            RAT=ODIF/T  
            QUAL=1./(1.+RAT**4) 
C        ...QW at obs locations:
            QOB(JOB)=WOB(JOB)*QUAL  
C        ...residuals (O-A) at obs locations:
            AO(JOB)=ODIF*QOB(JOB) 
  251    CONTINUE  
C     ...Spread QW and QW(O-A) to grid points:
         CALL OTOG2(XOB,YOB,NOB,NX,NY,QOB,AO,W,AI) 
C     ...Apply filter to gridded QW and QW(O-A) fields:
         CALL RF2D2(NX,NY,NSI,B1,B2,W,AI)
C     ...Update analysis according to modified eq(8):
         WMZA=WBG*WMA(IAI) 
C     ...Loop over grid points:
         DO 255 LXY=1,NXY  
            AA(LXY)=AA(LXY)+(AI(LXY)+WMZA*(AM(LXY)-AA(LXY)))
     *           /(WMZA+W(LXY))
C     ...Endloop over grid points:
  255    CONTINUE 
  300 CONTINUE  
      RETURN  
      END  

      SUBROUTINE RF2D2(NX,NY,NSI,B1,B2,A,B) 
c***********************************************************************
c* Function Name: rf2d2
c* Purpose: Filter
c* Usage: call RF2D2(NX,NY,NSI,B1,B2,A,B) 
c* Description: APPLY RECURSIVE FILTERS IN X AND Y ON A BOUNDED DOMAIN  
c*              TO 2 FIELDS, A AND B.  
c* Inputs:
c* Var_name                Description
c* --------                -----------
c* NX,NY                   No. of points in x,y
c* NSI                     No. of passes of filter
c* B1(NX,NY), B2(NX,NY)    Smoothing parameters b=1-alfa
c*                         for x, y filter
c* A(NX,NY), B(NX,NY)      Fields to be smoothed
c*
c* Outputs:
c* Var_name                Description
c* --------                -----------
c* A(NX,NY), B(NX,NY)      Fields to be smoothed
c***********************************************************************
      DIMENSION A(NX*NY),B(NX*NY),B1(NX*NY),B2(NX*NY) 
C     BEFORE WE BEGIN, A BRIEF NOTE ABOUT THE INDEXING CONVENTION: 
C     THE 2-D ARRAYS ARE MANIPULATED HERE AS 1-D ARRAYS FOR SPEEDINESS  
C     THE INDEX REPRESENTING GRID-POINT (IX,IY) IS, 
C     LXY = IX + (IY-1)*NX  
C     I.E., LXY = L0 + IX + LY*KY  
C     THE PART OF THIS EXPRESSION COMMON TO A LINE, E.G., Y  CONST, 
C     IS THE "LINE-INDEX", LY = L0 + LY*KY  
C     KY AND L0 ARE, OF COURSE, DEFINED BY: 
      KY=NX  
      L0=-KY  
      NX1=NX-1  
      NY1=NY-1 
C  ...Loop over no of passes:
      DO 200 ISI=1,NSI  
C     ...Loop over gridpoints in y:
         DO 201 IY=1,NY  
            LY=L0+IY*KY  
            LXY=1+LY  
            LXYP=LXY+1
C     SCAN FORWARDS AND BACKWARDS WITH THE BASIC X-FILTER  
C        ...Forward in x:
C        ...B.C. for i=1:
            BE=B1(LXY) 
            AL=1.-BE  
            AL2=AL*AL  
            AL3=AL*AL2  
            GA=1.-AL2  
            BEOGA=BE/GA  
            BEOGAS=BEOGA/GA  
            IF(ISI.EQ.1)THEN  
C     ...Eq(A1):
               A(LXY)=BE*A(LXY) 
               B(LXY)=BE*B(LXY) 
            ELSEIF(ISI.EQ.2)THEN  
C     ...Eq(A6):
               A(LXY)=BEOGA*A(LXY) 
               B(LXY)=BEOGA*B(LXY) 
            ELSEIF(ISI.GE.3)THEN  
C     ...Eq(A7) (mislabeled A6 in paper):
               A(LXY)=BEOGAS*(A(LXY)-AL3*A(LXYP))
               B(LXY)=BEOGAS*(B(LXY)-AL3*B(LXYP))
            ENDIF  
            DO 202 IX=2,NX  
               LXY=IX+LY  
               LXYM=LXY-1  
               BE=B1(LXY) 
               AL=1.-BE  
C           ...Eq(1):
               A(LXY)=BE*A(LXY)+AL*A(LXYM) 
               B(LXY)=BE*B(LXY)+AL*B(LXYM) 
  202       CONTINUE
C        ...Backward in x:
C        ...B.C. for i=NX:
            LXY=NX+LY  
            LXYM=LXY-1  
            BE=B1(LXY) 
            AL=1.-BE  
            AL2=AL*AL  
            AL3=AL*AL2  
            GA=1.-AL2  
            BEOGA=BE/GA  
            BEOGAS=BEOGA/GA  
            IF(ISI.EQ.1)THEN  
               A(LXY)=BEOGA*A(LXY) 
               B(LXY)=BEOGA*B(LXY)  
            ELSEIF(ISI.GE.2)THEN  
               A(LXY)=BEOGAS*(A(LXY)-AL3*A(LXYM))
               B(LXY)=BEOGAS*(B(LXY)-AL3*B(LXYM))
            ENDIF  
            DO 203 IX=NX1,1,-1  
               LXY=IX+LY  
               LXYP=LXY+1  
               BE=B1(LXY) 
               AL=1.-BE  
               A(LXY)=BE*A(LXY)+AL*A(LXYP) 
               B(LXY)=BE*B(LXY)+AL*B(LXYP) 
 203        CONTINUE  
C     ...Endloop over gridpoints in y:
  201    CONTINUE  
C     SCAN FORWARDS AND BACKWARDS WITH THE BASIC Y-FILTER  
C     ...Loop over gridpoints in x:
         DO 204 IX=1,NX
C        ...Forward in y:
C        ...B.C. for j=1:
            LX=IX+L0  
            LXY=KY+LX  
            LXYP=LXY+KY  
            BE=B2(LXY) 
            AL=1.-BE  
            AL2=AL*AL  
            AL3=AL*AL2  
            GA=1.-AL2  
            BEOGA=BE/GA  
            BEOGAS=BEOGA/GA  
            IF(ISI.EQ.1)THEN  
               A(LXY)=BE*A(LXY) 
               B(LXY)=BE*B(LXY) 
            ELSEIF(ISI.EQ.2)THEN  
               A(LXY)=BEOGA*A(LXY) 
               B(LXY)=BEOGA*B(LXY) 
            ELSEIF(ISI.GE.3)THEN  
               A(LXY)=BEOGAS*(A(LXY)-AL3*A(LXYP))
               B(LXY)=BEOGAS*(B(LXY)-AL3*B(LXYP))
            ENDIF  
            DO 205 IY=2,NY  
               LXY=IY*KY+LX  
               LXYM=LXY-KY  
               BE=B2(LXY) 
               AL=1.-BE  
               A(LXY)=BE*A(LXY)+AL*A(LXYM) 
               B(LXY)=BE*B(LXY)+AL*B(LXYM) 
  205       CONTINUE
C        ...Backward in y:
C        ...B.C. for j=NY:
            LXY=NY*KY+LX  
            LXYM=LXY-KY  
            BE=B2(LXY) 
            AL=1.-BE  
            AL2=AL*AL  
            AL3=AL*AL2  
            GA=1.-AL2  
            BEOGA=BE/GA  
            BEOGAS=BEOGA/GA  
            IF(ISI.EQ.1)THEN  
               A(LXY)=BEOGA*A(LXY) 
               B(LXY)=BEOGA*B(LXY) 
            ELSEIF(ISI.GE.2)THEN  
               A(LXY)=BEOGAS*(A(LXY)-AL3*A(LXYM))
               B(LXY)=BEOGAS*(B(LXY)-AL3*B(LXYM))
            ENDIF
            DO 207 IY=NY1,1,-1  
               LXY=LX+IY*KY  
               LXYP=LXY+KY  
               BE=B2(LXY) 
               AL=1.-BE  
               A(LXY)=BE*A(LXY)+AL*A(LXYP) 
               B(LXY)=BE*B(LXY)+AL*B(LXYP)
  207       CONTINUE  
C     ...Endloop over gridpoints in y:
  204    CONTINUE  
C  ...Endloop over passes:
  200 CONTINUE  
      RETURN  
      END  

      SUBROUTINE OTOG2(XOB,YOB,NOB,NX,NY,AOB,BOB,AF,BF) 
c***********************************************************************
c* Function Name: otog2
c* Purpose: Disperse observation at the observation points.
c* Usage: call OTOG2(XOB,YOB,NOB,NX,NY,AOB,BOB,AF,BF) 
c* Description: DISPERSE QOB AT OB POINTS TO QF AT GRID POINTS  
c*              THIS OPERATION IS THE ADJOINT OF THAT PERFORMED 
c*              IN SUBR. GTOO  
c* Inputs:
c* Var_name              Description
c* --------              -----------
c* XOB(NOB), YOB(NOB)    Normalized x, y coordinates of
c*                       obs (see comment in RFAN1Z)
c* NOB                   No. of observations
c* NX, NY                No. of gridpoints in x, y
c* AOB(NOB), BOB(NOB)    Field 1 and 2 at obs locations
c*
c* Outputs:
c* Var_name              Description
c* --------              -----------
c* AF(NX,NY), BF(NX,NY)  Field 1 and 2 at gridpoints
c***********************************************************************
C  ...All 2-D arrays are treated as 1-D arrays locally, and array 
C  ...indices computed accordingly
      DIMENSION XOB(NOB),YOB(NOB),AOB(NOB),BOB(NOB),AF(NX*NY),BF(NX*NY) 
      KY=NX  
      L0=-KY  
C  ...Loop over observations:
      DO 200 IOB=1,NOB  
         if(xob(iob).le.0) goto 200
C     ...Find surrounding 4 gridpoints: 
C        L11=(IX1,IY1)    L21=(IX1+1,IY1)
C        L12=(IX1,IY1+1)  L22=(IX1+1,IY1+1)
         IX=XOB(IOB) 
         RX=XOB(IOB)-IX  
         IX1=IX+1  
         IY=YOB(IOB) 
         IY=MIN(IY,NY-2)
         RY=YOB(IOB)-IY  

         IY1=IY+1  
         L11=L0+IX1+IY1*KY  
         L21=L11+1  
         L12=L11+KY  
         L22=L12+1  
C     ...Apply Equation shown in section 2.2:
         RX1=1.-RX  
         RY1=1.-RY  
         RX1RY1=RX1*RY1  
         RXRY1=RX*RY1  
         RX1RY=RX1*RY  
         RXRY=RX*RY  
         AIOB=AOB(IOB) 
         BIOB=BOB(IOB) 
C     ...Field 1:
         AF(L11)=AF(L11)+RX1RY1*AIOB  
         AF(L12)=AF(L12)+RX1RY*AIOB  
         AF(L21)=AF(L21)+RXRY1*AIOB  
         AF(L22)=AF(L22)+RXRY*AIOB  
C     ...Field 2:
         BF(L11)=BF(L11)+RX1RY1*BIOB  
         BF(L12)=BF(L12)+RX1RY*BIOB  
         BF(L21)=BF(L21)+RXRY1*BIOB  
         BF(L22)=BF(L22)+RXRY*BIOB  
C  ...Endloop over observations:
  200 CONTINUE  
      RETURN  
      END  

      SUBROUTINE GTOO(XOB,YOB,NOB,NX,NY,QF,QOB) 
c***********************************************************************
c* Function Name: gtoo
c* Purpose: Interpolation
c* Usage: call GTOO(XOB,YOB,NOB,NX,NY,QF,QOB)
c* Description: INTERPOLATE QF AT GRID POINTS TO QOB AT OB POINTS  
c*              Algorithm- Bilinear interpolation
c* 
c* Inputs:
c* Var_name              Description
c* --------              -----------
c* XOB(NOB), YOB(NOB)    Normalized x, y coordinates of
c*                       obs (see comment in RFAN1Z)
c* NOB                   No. of observations
c* NX, NY                No. of gridpoints in x, y
c* QF(NX,NY)             Field at gridpoints
c*
c* Outputs:
c* Var_name              Description
c* --------              -----------
c* QOB(NOB)              Field at obs locations
c***********************************************************************
C  ...All 2-D arrays are treated as 1-D arrays locally, and array 
C  ...indices computed accordingly
      DIMENSION QOB(NOB),XOB(NOB),YOB(NOB),QF(NX*NY) 
      KY=NX  
      L0=-KY  
C  ...Loop over observations:
      DO 200 IOB=1,NOB 
         if(xob(iob).le.0.) goto 200
C     ...Find surrounding 4 gridpoints: 
C        L11=(IX1,IY1)    L21=(IX1+1,IY1)
C        L12=(IX1,IY1+1)  L22=(IX1+1,IY1+1)
         IX=XOB(IOB) 
         RX=XOB(IOB)-IX  
         IX1=IX+1  
         IY=YOB(IOB) 
         IY=MIN(IY,NY-2)
         RY=YOB(IOB)-IY  
         IY1=IY+1  
         L11=L0+IX1+IY1*KY  
         L21=L11+1  
         L12=L11+KY  
         L22=L12+1
         RX1=1.-RX  
         RY1=1.-RY 
C     ...Apply bilinear interpolation formula:
         RX1RY1=RX1*RY1  
         RXRY1=RX*RY1  
         RX1RY=RX1*RY  
         RXRY=RX*RY  
         QOB(IOB)=QF(L11)*RX1RY1  
     *        +QF(L12)*RX1RY  
     *        +QF(L21)*RXRY1  
     *        +QF(L22)*RXRY
C  ...Endloop over observations:
  200 CONTINUE  
      RETURN  
      END
 
