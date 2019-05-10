c-------------------------------------------------------------------------------
c     $Name$ 
c     $Id$ 
c     Copyright AER, Inc., 2002, 2003. All rights Reserved.
c-------------------------------------------------------------------------------

      subroutine fft (z,mval,iwk)                                       
************************************************************************
c     E2.4                                                              
c Uses double precision  Rev Date May 9,1987   
c                                                                       
c   FUNCTION            - compute the fast fourier transform, given a   
c                           complex vector of length equal to a power   
c                           of two                                      
c   USAGE               - call 
c   PARAMETERS  z       - complex vector of length n=2**m               
c                           which contains on input the                 
c                           data to be transformed. on                  
c                           output,a contains the fourier               
c                           coe
c                m      - n = 2**m is the number of data points.        
c                         m= +n 
c                             transform.                                
c                         m= -n 
c                             transform.                                
c                iwk    - work area vector of length m+1.               
c   PRECISION           - single                                        
c   LANGUAGE            - fortran                                       
c   LATEST REVISION     - April 16, 1980                                
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension iwk(*),z(*),z0(2),z1(2),z2(2),z3(2)                     
      complex*16 z,za0,za1,za2,za3,ak2                                  
      equivalence (za0,z0(1)),(za1,z1(1)),(za2,z2(1)),(za3,z3(1)),      
     &     (a0,z0(1)),(b0,z0(2)),(a1,z1(1)),(b1,z1(2)),                 
     &     (a2,z2(1)),(b2,z2(2)),(a3,z3(1)),(b3,z3(2))                  
      data               sq,sk,ck /.70710678118655d0,.38268343236509d0, 
     &                            .92387953251129d0/                    
      data               twopi/6.2831853071796d0/,zero/0.0d0/,one/1.0d0/
c     sq=sqrt2/2,sk=sin(pi/8),ck=cos(pi/8) 
c     twopi=2*pi                           
c     the variable sign (also a fortran 90 intinsic) is replaced with ssign
      ssign = 1.0d0                                                      
      if (mval .lt. 0) ssign = -1.0d0                                    
      m= iabs(mval)                                                     
      mp = m+1                                                          
      n = 2**m                                                          
      if (ssign .eq. 1.0d0) go to 4                                      
      do 2 i=1,n                                                        
         z(i) = conjg(z(i))                                            
 2    continue                                                          
 4    iwk(1)=1                                                          
      mm = (m/2)*2                                                      
      kn = n+1                                                          
c-- initialize work vector               
      do 5 i=2,mp                                                       
         iwk(i) = iwk(i-1)+iwk(i-1)                                     
 5    continue                                                          
      rad=twopi/n                                                       
      mk = m - 4                                                        
      kb = 1                                                            
      if (mm .eq. m) go to 15                                           
      k2 = kn                                                           
      k0 = iwk(mm+1) + kb                                               
 10   k2 = k2 - 1                                                       
      k0 = k0 - 1                                                       
      ak2 = z(k2)                                                       
      z(k2) = z(k0)- ak2                                                
      z(k0) = z(k0)+ ak2                                                
      if (k0 .gt. kb) go to 10                                          
 15   c1 = one                                                          
      s1 = zero                                                         
      jj = 0                                                            
      k = mm - 1                                                        
      j = 4                                                             
      if (k .ge. 1) go to 30                                            
      go to 9005                                                        
 20   if (iwk(j) .gt. jj) go to 25                                      
      jj = jj - iwk(j)                                                  
      j = j-1                                                           
      if (iwk(j) .gt. jj) go to 25                                      
      jj = jj - iwk(j)                                                  
      j = j - 1                                                         
      k = k + 2                                                         
      go to 20                                                          
 25   jj = iwk(j) + jj                                                  
      j = 4                                                             
 30   isp = iwk(k)                                                      
      if (jj .eq. 0) go to 40                                           
c-- reset trigonometric parameters       
      c2 = jj * isp * rad                                               
      c1 = cos(c2)                                                      
      s1 = sin(c2)                                                      
 35   c2 = c1 * c1 - s1 * s1                                            
      s2 = c1 * (s1 + s1)                                               
      c3 = c2 * c1 - s2 * s1                                            
      s3 = c2 * s1 + s2 * c1                                            
 40   jsp = isp + kb                                                    
c-- determine fourier coef. in groups of 4                     
      do 50 i=1,isp                                                     
         k0 = jsp - i                                                   
         k1 = k0 + isp                                                  
         k2 = k1 + isp                                                  
         k3 = k2 + isp                                                  
         za0 = z(k0)                                                    
         za1 = z(k1)                                                    
         za2 = z(k2)                                                    
         za3 = z(k3)                                                    
         if (s1 .eq. zero) go to 45                                     
         temp = a1                                                      
         a1 = a1 * c1 - b1 * s1                                         
         b1 = temp * s1 + b1 * c1                                       
         temp = a2                                                      
         a2 = a2 * c2 - b2 * s2                                        
         b2 = temp * s2 + b2 * c2                                       
         temp = a3                                                      
         a3 = a3 * c3 - b3 * s3                                         
         b3 = temp * s3 + b3 * c3                                       
 45      temp = a0 + a2                                                 
         a2 = a0 - a2                                                   
         a0 = temp                                                      
         temp = a1 + a3                                                 
         a3 = a1 - a3                                                   
         a1 = temp                                                      
         temp = b0 + b2                                                 
         b2 = b0 - b2                                                   
         b0 = temp                                                      
         temp = b1 + b3                                                 
         b3 = b1 - b3                                                   
         b1 = temp                                                      
         z(k0) = dcmplx(a0+a1,b0+b1)                                    
         z(k1) = dcmplx(a0-a1,b0-b1)                                    
         z(k2) = dcmplx(a2-b3,b2+a3)                                    
         z(k3) = dcmplx(a2+b3,b2-a3)                                    
 50   continue                                                          
      if (k .le. 1) go to 55                                            
      k = k - 2                                                         
      go to 30                                                          
 55   kb = k3 + isp                                                     
c-- check for completion of final iteration                          
      if (kn .le. kb) go to 9005                                        
      if (j .ne. 1) go to 60                                            
      k = 3                                                             
      j = mk                                                           
      go to 20                                                          
 60   j = j - 1                                                         
      c2 = c1                                                           
      if (j .ne. 2) go to 65                                            
      c1 = c1 * ck + s1 * sk                                            
      s1 = s1 * ck - c2 * sk                                            
      go to 35                                                          
 65   c1 = (c1 - s1) * sq                                               
      s1 = (c2 + s1) * sq                                               
      go to 35                                                          
 9005 continue                                                          
      if (ssign .eq. 1.0) go to 75                                       
      xn = n                                                            
      do 70 i=1,n                                                       
      z(i)=conjg(z(i))/xn                                              
 70   continue                                                          
 75   call qxz136(z,m,iwk)                                              
      return                                                           
      end                                                             
c --- subprogram qxz136 --- formerly known as routine  ffrdr2  ---   
      subroutine qxz136 (z,m,iwk)                                   
c-ffrdr2--------s-------library 3---------------------------------------
c                                                                       
************************************************************************
c                                                                       
c   FUNCTION            - This subroutine permutes a complex data vector
c                           in reverse binary order to normal order. The
c                           routine can also be used to permute a com-  
c                           plex data vector in normal order to reverse 
c                           binary order since the permutation is sym-  
c                           metric.                                     
c   USAGE               - call qxz136(z,m,iwk)         
c   PARAMETERS  z       - complex vector of length n=2**m which         
c                           contains on input the data to be            
c                           permuted. on output, z contains the         
c                           permuted data vector.                       
c                m      - n=2**m is the number of data points.          
c                iwk    - work area vector of length m+1                
c   PRECISION           - single                                        
c   LANGUAGE            - fortran                                       
c   LATEST REVISION     - March 16, 1973                                
c     ----------------------------------------------------------------- 
c                                                                       
      implicit real*8 (a-h,o-z)
      dimension iwk(*),z(*)                                             
      complex*16 z,temp                                                 
c                                                                       
      if(m .le. 1) go to 9005                                           
      mp = m+1                                                          
      jj = 1                                                            
c-- initialize work vector               
      iwk(1) = 1                                                        
      do 5  i = 2,mp                                                    
         iwk(i) = iwk(i-1) * 2                                          
 5    continue                                                          
      n4 = iwk(mp-2)                                                    
      if (m .gt. 2) n8 = iwk(mp-3)                                     
      n2 = iwk(mp-1)                                                    
      lm = n2                                                           
      nn = iwk(mp)+1                                                    
      mp = mp-4                                                         
c-- determine indices and switch a*s     
      j = 2                                                             
 10   jk = jj + n2                                                      
      temp= z(j)                                                        
      z(j)=z(jk)                                                        
      z(jk)= temp                                                       
      j = j+1                                                           
      if (jj .gt. n4) go to 15                                          
      jj = jj + n4                                                      
      go to 35                                                          
 15   jj = jj - n4                                                      
      if (jj .gt. n8) go to 20                                          
      jj = jj + n8                                                      
      go to 35                                                          
 20   jj = jj - n8                                                      
      k = mp                                                            
 25   if (iwk(k) .ge. jj) go to 30                                      
      jj = jj - iwk(k)                                                  
      k = k - 1                                                         
      go to 25                                                          
 30   jj = iwk(k) + jj                                                  
 35   if (jj .le. j) go to 40                                           
      k = nn - j                                                        
      jk = nn - jj                                                      
      temp= z(j)                                                        
      z(j) = z(jj)                                                      
      z(jj) = temp                                                      
      temp = z(k)                                                       
      z(k) = z(jk)                                                      
      z(jk)= temp                                                       
 40   j = j + 1                                                         
c-- cycle repeated until limiting number 
c-- of changes is achieved             
      if (j .le. lm) go to 10                                           
 9005 return                                                            
      end                                                              

