c     c********************************************************************
      subroutine prof2CLW(P,PROF,NLEV,PSFC,QINT,CTOP,CTHK)

C     CALCULATES THE TOTAL INTEGRATED CLOUD LIQUID WATER IN A 
C     COLUMN IN KG PER SQUARE METER, PLUS CLOUD TOP PRESSURE AND
C     CLOUD THICKNESS PRESSURE.
c     JFG:  Converted to use clw in kg/kg (was g/kg)

      DIMENSION P(*),PROF(*)

C     LOCATE SURFACE
C     Always used 1000 mb Q for Qsfc since the Qsfc for land
C     points, prof(iqsfc), is not correct.

      QSFC=PROF(1)   
      DO I=2,nlev
         IF(PSFC.GT.P(I)) exit
      END DO
C     INITIALIZE
      QINT=0.
      NN=I
      P0=PSFC
      Q0=QSFC

c     estimate qsfc when Psfc < 1000mb

C     For clw need to do this interpolation linearly
c     if(nn.gt.2)then
c     c         x=alog(prof(nn-1+15)/prof(nn-2+15))/alog(p(nn-1)/p(nn-2))
c     c         q0=prof(nn-2+15)*(psfc/p(nn-2))**x
c     x=alog(prof(nn)/prof(nn-1))/alog(p(nn)/p(nn-1))
c     q0=prof(nn-1)*(psfc/p(nn-1))**x
c     print*,prof(nn),prof(nn-1)
c     print *,'tiwv: ',psfc,q0,prof(1)
c     end if

c     Only go up to 100 mb.
      do i=1,nlev
         if (p(i) <= 100.0) then
            imax = i
            exit
         endif
      enddo

      nbase=0
      DO I=NN,imax                ! only go up to 100 mb.
         P1=P(I)
         Q1=PROF(I)
         Q1=MAX(0.0,Q1)
c     JFG:  Added E-3 factor (and below) to account for kg/kg clw
         if(q1.gt.0.01E-3.and.nbase.eq.0) nbase=i
C     X=ALOG(Q1/Q0)/ALOG(P1/P0)
C     QINT=QINT+1/(9.80*(X+1))*(P0*Q0-P1*Q1)/10.
CJFG was:
C         QINT=QINT+.5*(Q1+Q0)*(P0-P1)/98.
CJFG now is:
         QINT=QINT+.5*(Q1+Q0)*(P0-P1)*100.0/9.8
         P0=P1
         Q0=Q1
      END DO

      ntop=0
      do i=imax,nn,-1		! only go up to 100 mb
         if(prof(i).gt.0.01E-3) then
            ntop=i
            exit
         end if
      end do
c      ctop=p(ntop)
c      cthk=p(nbase) - ctop

c      if(qint.lt.0.01) then
c     JFG: Changed to be consistent with base/top criteria
c     And, yes, ctop and cthk are in pressure units
      if (ntop.ne.0.and.nbase.ne.0) then
         if (ntop.eq.nbase) then
            ntop = nbase+1
         endif
         ctop = p(ntop)
         cthk = p(nbase) - ctop
      else
         ctop=500.0
         cthk=50.0
         qint=0.0
      end if
C     
 100  RETURN
      END
