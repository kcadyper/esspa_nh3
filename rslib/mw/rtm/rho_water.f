C-----------------------------------------------------------------------C
      function RHOW(TC)
C     density of water in g/cm**3
C     input parameter TC is temperature in deg. celcius
      integer I
      real RHOW,TC,A(0:5),B,TEMP
      data A/999.8396,18.224944,-7.92221e-3,-55.44846e-6,
     &        149.7562e-9,-393.2952e-12/
      data B/18.159725e-3/
      TEMP=1.0
      RHOW=A(0)
      do I=1,5
        TEMP = TEMP*TC
        RHOW=RHOW + A(I)*TEMP
      end do
      RHOW = RHOW*(1.0e-3)/(1.0+B*TC)
      return
      end
C-----------------------------------------------------------------------C
      function RHOI(TC)
C     density of ice in g/cm**3
C     input parameter TC is temperature in deg. celsius
      real A0,A1,A2,TC,RHOI
      data A0,A1,A2/0.9167,-1.75e-4,-5.0e-7/
      RHOI=A0+A1*TC+A2*TC*TC
      return
      end
C-----------------------------------------------------------------------C

