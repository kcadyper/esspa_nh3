      subroutine DIE2NDX(ER,EI,NR,NI)
C     real and imaginary compenents of the dielectric constant (ER,EI)
C     are converted into the real and imaginary compenents of the
C     index of refraction.   The equation used is 9.2 from Bohren and 
C     Huffman "Absorption and Scattering of Light by Small Particles"
C     (1983)
C     July 1998 : B.T. Johnson - Purdue University - Department of 
C     Earth and Atmospheric Sciences.  jbenjam@meso.eas.purdue.edu
      implicit none
c      real NR,NI,ER,EI
c      complex E,N

c     **** 1 Jul 04 test:
      real NI,ER,EI
      complex E,N,NR
c     ******************

       E=CMPLX(ER,EI)
       N=CSQRT(E)
c       NR=SNGL(N)
       NR=REAL(N)
       NI=AIMAG(N)

      return
      end
