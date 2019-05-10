      subroutine NDX2DIE(NR,NI,ER,EI)
C************************************************************************
C
C     Purpose : Convert index of refraction to dielectric constant.
C 
C     Real and imaginary components of the index of refraction (NR,NI)
C     are converted into the real and imaginary components
C     of the dielectric function.  The equation used is 9.1 from Bohren
C     and Huffman "Absorption and Scattering of Light by Small Particles"
C     (1983)
C
C     Inputs   NR = real part of index of refraction
C              NI = imaginary part of index of refraction
C
C     Outputs  ER = real part of dielectric function
C              EI = imaginary part of dielectric function
C
C     July 2000 : Benjamin T. Johnson - University of Wisconsin
C                 [jbenjam@aos.wisc.edu]
C*********************************************************************

      implicit none
      real NR,NI,ER,EI

      ER = NR**2-NI**2
      EI = 2*NR*NI

      return
      end 
