      program test_gauss
      implicit none

      integer, parameter :: nf=100000
      real foo(nf)
      integer(kind=4) i, iseed, iseed_start
      real gauss
      real s, am
      real ave, adev, sdev, var, skew, curt
      iseed_start = 123456
      iseed_start= 0

      iseed = iseed_start

      s = 2.5
      am = 1.3
      print *,'iseed, nf, s, am: ', iseed_start, nf,s, am
      do i=1,nf
         foo(i)=gauss(iseed, s, am)
         iseed=0
      end do

      call moment(foo, nf, ave, adev, sdev, var, skew, curt)

      print *,'ave:  ', ave
      print *,'adev: ', adev
      print *,'sdev: ', sdev
      print *,'var:  ', var
      print *,'skew: ', skew
      print *,'curt: ', curt

      open(10, file='random_numbers.dat', status='unknown',
     &     form='formatted')

      write(10,*) nf, iseed_start, s,am
c*****Need to format record, otherwise, on sun,  it will write all nf numbers
c     to a single record
      write(10,fmt='(5g15.7)') foo
      close (10)

      end

      SUBROUTINE MOMENT(DATA,N,AVE,ADEV,SDEV,VAR,SKEW,CURT)
      DIMENSION DATA(N)
      IF(N.LE.1)PAUSE 'N must be at least 2'
      S=0.
      DO 11 J=1,N
        S=S+DATA(J)
11    CONTINUE
      AVE=S/N
      ADEV=0.
      VAR=0.
      SKEW=0.
      CURT=0. 
      DO 12 J=1,N
        S=DATA(J)-AVE
        ADEV=ADEV+ABS(S)
        P=S*S
        VAR=VAR+P
        P=P*S
        SKEW=SKEW+P
        P=P*S
        CURT=CURT+P
12    CONTINUE
      ADEV=ADEV/N
      VAR=VAR/(N-1)
      SDEV=SQRT(VAR)
      IF(VAR.NE.0.)THEN
        SKEW=SKEW/(N*SDEV**3)
        CURT=CURT/(N*VAR**2)-3.  
      ELSE
        print *, 'no skew or kurtosis when zero variance'
      ENDIF
      RETURN
      END
