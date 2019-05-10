      subroutine compress(dcov,nparg,ntran,u_mw,ucov)

! Use dynamic memory allocation for local variables.


! NOTE:  When the covariance program is used to generate a covariance for
! all parameters may want to call calc_eigenvector once each for temperature
! water vapor and emissivity.  The break points in the the background vector
! and covariance matrix could be passed in to compress ( even the number of
! break points to make it completely general) and just the sub dimensions
! of each array could be passed into calc_eigenvector.

      implicit none

      integer AllocateStatus
      integer :: nparg,ntran
      real  dcov(nparg,nparg),ucov(ntran,ntran)
      real u_mw(nparg,ntran)
      integer :: i

!-- local vars
!      real tmp1(nparg*nparg)

! declare allocatable for dynamic or run-time memorary allocation
      real, dimension(:,:), allocatable :: tmp1

      allocate (tmp1(ntran,nparg),stat=AllocateStatus) 
      if (AllocateStatus /=0)                                           &
     &     STOP "** ERROR in compress: Not enough memory tmp1 **"

!      do k=1,nparg
!         print*,k,dcov(k,k)
!      end do

!---  eigenvector transformation of mw covariance matrix:
      
      tmp1=matmul(transpose(u_mw),dcov)
      
      ucov=matmul(tmp1,u_mw)

!      do i =1, ntran
!        print *,ucov(i,i)
!      enddo

!--- Need to use something like map_geo2retr and dback, un, and u_mw to get umean or
!--- background in retrieval space.  Will insert here.  Need to determine exactly what is
!--- is expected by retrieval code.

!--- Note any correlation between parameters such as t and w should be 
!--- retained in dcov which influences ucov through the product of ut_mw and dcov.
!--- Question: does ut_mw have all zeros on the diagonals?  Ask Xu.


      deallocate(tmp1)
      
      return
      end

!--------------------------------------------------------------------------
!-- compute the eigenvectors for a block along the diagonal of a covariance
!--------------------------------------------------------------------------
      subroutine calceigenvector(dcov,u_mw,wr,nparg,np1,np2)

!***********************************************************************
!* Function name: calceigenvector
!* Description: calculate eigenvectors.
!* Inputs:
!* Var_Name      I*n          Description
!* --------      ---          -----------
!* dcov            real         covariance matrix.
!*
!* Outputs:
!* Var_Name      I*n        Description
!* --------      ---        -----------
!* u_mw          real       eigenvector matrix.
!* wr            real       eigenvalue vector
!*
!*  Developed by Atmospheric and Environmental Research, Inc.            
!*  Copyright: Atmospheric and Environmental Research, Inc., 1999        
!***********************************************************************

      implicit none

! Internal math is double precision
      integer, parameter :: dbl=selected_real_kind(12)
      
!-- Input
      real dcov(nparg,nparg)
      integer :: nparg,np1,np2

!-- Output

      real u_mw(nparg,nparg),wr(nparg)

!-- Local

      integer AllocateStatus
!      real w(npar1),usem(nparg,nparg),vsem(nparg,nparg)
!      real vnem(nparg,nparg)

!      real, dimension(:), allocatable :: w
!      real, dimension(:,:), allocatable :: usem,vsem,vnem,unem
      real(kind=dbl), dimension(:), allocatable :: w
      real(kind=dbl), dimension(:,:), allocatable :: usem,vsem,vnem,unem
      integer :: npdim, npdim1
      integer :: i,j,ki,kj,itmp
      real(kind=dbl) :: tmp

      npdim=np2-np1+1

      npdim1=npdim+1

!      do i=1,nparg
!         print*,i,dcov(i,i)
!      end do

      allocate (w(npdim1),stat=AllocateStatus) 
      if(AllocateStatus /=0) STOP "** Not enough memory w **"

      allocate (usem(npdim,npdim),vsem(npdim,npdim),                    &
     &     vnem(npdim,npdim),unem(npdim,npdim),                         &
     &     stat=AllocateStatus)
 
      if(AllocateStatus /=0) STOP "** Not enough memory usem **"

!---  eigenvector cov conversion:

!------get eigenvectors and assign them to u_mw---

      do i=1,npdim
         do j=1,npdim
            usem(i,j)=0.
            vsem(i,j)=0.
            unem(i,j)=0.0
            vnem(i,j)=0.0
         end do
      end do

      w(1:npdim1)=0.0

      ki=0
      do i=np1,np2
         ki=ki+1
         kj=0
         do j=np1,np2
            kj=kj+1
            usem(ki,kj)=real(dcov(i,j),kind=dbl)
         enddo
      enddo   


      call svdcmp(usem,npdim,npdim,npdim,npdim,w,vsem)


!      print *,'-----'
!      print *, w
!      print *,'-----'

      call svdsort(usem,vsem,w,unem,vnem,npdim,npdim,npdim,npdim)
       
!      print *,'-----'
!      print *, w
!      print *,'-----'

      ki=0
      do i=np1,np2
         ki=ki+1
         kj=0
         do j=np1,np2
            kj=kj+1
            u_mw(i,j) = real(unem(ki,kj))
         end do
         wr(i)=w(ki)
      end do

      deallocate(w)
      deallocate(usem,vsem,vnem,unem)

      return
      end
  


!******************** SVD Routines ******************************



      subroutine svdcmp(a,m,n,mp,np,w,v)
!************************************************************************
!* Function Name: svdcmp                                                  
!* Description: Algorithm for constructing the singular value           
!*              decomposition of any matrix.                            
!*              Given a matrix a(1:m,1:n), with phyiscal dimensions     
!*              mp by np, this routine computes its singular value      
!*              decomposition, A=U*W*V^T. The matrix U replaces a       
!*              on output. The diagonal matrix of singular values       
!*              W is output as a vector w(1:n). The matrix V (not       
!*              the transpose V^T) is output as v(1:n,1:n).             
!* Usage: call svdcmp(a,m,n,mp,np,w,v)
!* Input Args:    I*n    Var_Name     Description                       
!*                real     a          matrix                            
!*                integer  m,n        logical dimensions of a           
!*                integer  mp,np      physical dimensions of a          
!*                                                                       
!* Output Args:   I*n    Var_Name     Description                       
!*                real     w          1-by-n vector                     
!*                real     v          n-by-n square matrix              
!*                                                                       
!* Externals:     File_Name           Description                       
!*                pythag              function to calculate             
!*                                    sqrt(a**2+b**2)                 
!* Reference: Numerical Recipes in Fortran, Press et al, 1992           
!* Modified by Atmopsheric and Environmental Research, Inc., 1998       
!************************************************************************

      implicit none

      integer, parameter :: dbl=selected_real_kind(12)

      integer m,mp,n,np,nmax
      parameter (nmax=2500)      !Maximum anticipated value of n

      real(kind=dbl) a(mp,np),v(np,np),w(np+1)
      
!-- Uses pythag      
      integer i,its,j,jj,k,l,nm
      real(kind=dbl) anorm,c,f,g,h,s,scale,x,y,z,rv1(nmax),pythag

      g = 0.0
      scale = 0.0
!     z fix anorm from a amorn
      anorm = 0.0
      do 25 i = 1, n
         l = i + 1
         rv1(i) = scale*g
         g = 0.0
         s = 0.0
         scale = 0.0
         if (i.le.m) then
            do k = i, m
               scale = scale + abs(a(k,i))
            enddo
            if (scale.ne.0.0) then
               do k = i, m
                  a(k,i) = a(k,i)/scale
                  s = s + a(k,i)*a(k,i)
               enddo
               f = a(i,i)
               g = -sign(sqrt(s),f)
               h = f*g - s
               a(i,i) = f - g
               do j = l, n
                  s = 0.0
                  do k = i, m
                     s = s + a(k,i)*a(k,j)
                  enddo
                  f = s/h
                  do k = i, m
                     a(k,j) = a(k,j) + f*a(k,i)
                  enddo
               enddo
               do k =i, m
                  a(k,i) = scale*a(k,i)
               enddo
            endif
         endif
         w(i) = scale*g
         g = 0.0
         s = 0.0
         scale = 0.0
         if ((i.le.m).and.(i.ne.n)) then
            do k = l, n
               scale = scale + abs(a(i,k))
            enddo
            if (scale.ne.0.0) then
               do k = l, n
                  a(i,k) = a(i,k)/scale
                  s = s + a(i,k)*a(i,k)
               enddo
               f = a(i,l)
               g = -sign(sqrt(s),f)
               h = f*g - s
               a(i,l) = f - g
               do k = l, n
                  rv1(k) = a(i,k)/h
               enddo
               do j = l, m
                  s = 0.0
                  do k = l, n
                     s = s + a(j,k)*a(i,k)
                  enddo
                  do k = l, n
                     a(j,k) = a(j,k) + s*rv1(k)
                  enddo
               enddo
               do k = l, n
                  a(i,k) = scale*a(i,k)
               enddo
            endif
         endif
         anorm = max(anorm,(abs(w(i))+abs(rv1(i))))
  25     continue
      do 32 i = n, 1, -1 !Accumulation of right-hand transformation
         if (i.lt.n) then
            if (g.ne.0.0) then
               do j = l, n
                  v(j,i) = (a(i,j)/a(i,l))/g
               enddo
               do j = l, n
                  s = 0.0
                  do k = l, n
                     s = s + a(i,k)*v(k,j)
                  enddo
                  do k = l, n
                     v(k,j) = v(k,j) + s*v(k,i)
                  enddo
               enddo
            endif
            do j = l, n
               v(i,j) = 0.0
               v(j,i) = 0.0
            enddo
         endif
         v(i,i) = 1.0
         g = rv1(i)
         l = i
  32     continue
      do 39 i = min(m,n), 1, -1 !Accumulation of left-hand transformation
         l = i + 1
         g = w(i)
         do j = l, n
            a(i,j) = 0.0
         enddo
         if (g.ne.0.0) then
            g = 1.0/g
            do j = l, n
               s = 0.0
               do k = l, m
                  s = s + a(k,i)*a(k,j)
               enddo
               f = (s/a(i,i))*g
               do k = i, m
                  a(k,j) = a(k,j) + f*a(k,i)
               enddo
            enddo
            do j = i, m
               a(j,i) = a(j,i)*g
            enddo
         else
            do j = i, m
               a(j,i) = 0.0
            enddo
         endif
         a(i,i) = a(i,i) + 1.0
  39     continue
      do 49 k = n, 1, -1  !Diagonalization of bidiagonal form: Loop
         do its = 1, 100  !  over singular values, and over allowed 
                          !  iterations.
            do l = k, 1, -1     !Test for splitting.
               nm = l - 1       !Note that rv1(1) is always zero.

!     z fix abs from cab
!     z add new convergence test to get around not going to zero.
!     z inc the number of convergence iterations to 200

!     z if ((abs(rv1(l))+anorm).eq.anorm) goto 2               
               if (abs(rv1(l)).lt.1.e-20) goto 2
               if ((abs(w(nm))+anorm).eq.anorm) goto 1
            enddo
  1                   c = 0.0  !Cancellation of rv1(l), if l > 1
            s = 1.0
            do i = l, k
               f = s*rv1(i)
               rv1(i) = c*rv1(i)
               if ((abs(f)+anorm).eq.anorm) goto 2
               g = w(i)
               h = pythag(f,g)
               w(i) = h
               h = 1.0/h
               c = g*h
               s = -f*h
               do j = 1, m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = y*c + z*s
                  a(j,i) = -y*s + z*c
               enddo
            enddo
  2                   z = w(k)
            if (l.eq.k) then    !Convergence
               if (z.lt.0.0) then
                  w(k) = -z
                  do j = 1, n
                     v(j,k) = -v(j,k)
                  enddo
               endif
               goto 3
            endif
            if (its.eq.100) stop 'no convergence in svdcmp'
            x = w(l)            !Shift from bottom 2-by-2 minor.
            nm = k - 1
            y = w(nm)
            g = rv1(nm)
            h = rv1(k)
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g = pythag(f,1.0_dbl)
            f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
            c = 1.0             !Next QR transformation:
            s = 1.0
            do j = l, nm
               i = j + 1
               g = rv1(i)
               y = w(i)
               h = s*g
               g = c*g
               z = pythag(f,h)
               rv1(j) = z
               c = f/z
               s = h/z
               f = x*c + g*s
               g = -x*s + g*c
               h = y*s
               y = y*c
               do jj = 1, n
                  x = v(jj,j)
                  z = v(jj,i)
                  v(jj,j) = x*c + z*s
                  v(jj,i) = -x*s + z*c
               enddo
               z = pythag(f,h)
               w(j) = z         !Rotation can be arbitrary if z = 0.
               if (z.ne.0.0) then
                  z = 1.0/z
                  c = f*z
                  s = h*z
               endif
               f = c*g + s*y
               x = -s*g + c*y
               do jj = 1, m
                  y = a(jj,j)
                  z = a(jj,i)
                  a(jj,j) = y*c + z*s
                  a(jj,i) = -y*s + z*c
               enddo
            enddo
            rv1(l) = 0.0
            rv1(k) = f
            w(k) = x
         enddo
  3             continue
  49            continue

      return
      end

      function pythag(a,b)
!----------------------------------------------------------------------
!-- Computes sqrt(a**2+b**2) without destrutive underflow or overflow.
!----------------------------------------------------------------------

      implicit none

      integer, parameter :: dbl=selected_real_kind(12)

      real(kind=dbl) a,b
      real(kind=dbl) pythag
      real(kind=dbl) absa,absb

      absa = abs(a)
      absb = abs(b)
      if (absa.gt.absb) then
         pythag = absa*sqrt(1.0+(absb/absa)**2)
      else
         if (absb.eq.0.0) then
            pythag = 0.0
         else
            pythag = absb*sqrt(1.0+(absa/absb)**2)
         endif
      endif

      return
      end


      subroutine svdsort(u,v,w,unew,vnew,m,n,mp,np)
!************************************************************************
!* Function Name: svdsort
!* Description: sorting the matrics by descending order.
!* Usage: call svdsort(u,v,w,unew,vnew,m,n,mp,np)
!* Input Args:   
!* Var_Name     I*n          Description                       
!* --------     ---          ----------- 
!* u,v          real         matrices to be sorted.
!* w            real         sorting criterion.
!* m,n,mp,np    integer      dimensions.
!*
!* Output Args:   
!* Var_Name     I*n          Description                       
!* --------     ---          ----------- 
!* unew,vnew    real         matrices sorted.
!*                                                                       
!* Developed by Atmopsheric and Environmental Research, Inc., 1998       
!* CopyRight,  Atmopsheric and Environmental Research, Inc., 1998 
!************************************************************************

      implicit none

      integer, parameter :: dbl=selected_real_kind(12)

      integer m,mp,n,np,nmax
      parameter (nmax=2500)

      integer indx(nmax)
!      real u(mp,np),v(np,np),w(np)
!      real unew(mp,np),vnew(np,np)
      real(kind=dbl) u(mp,np),v(np,np),w(np+1)
      real(kind=dbl) unew(mp,np),vnew(np,np)
      integer :: i,j,itmp
      real(kind=dbl) :: tmp

      do i=1,n
         indx(i)=i
      enddo

      do i=1,n-1
         do j=1,n-i
            if(w(j) .lt. w(j+1)) then
               tmp=w(j)
               w(j)=w(j+1)
               w(j+1)=tmp
               itmp=indx(j)
               indx(j)=indx(j+1)
               indx(j+1)=itmp
            endif
         enddo
      enddo

      do i=1,n
         do j=1,m
            unew(j,i)=u(j,indx(i))
         enddo
         do j=1,n
            vnew(j,i)=v(j,indx(i))
         enddo
      enddo

      return
      end


      subroutine svbksb(u,w,v,m,n,mp,np,b,x)
!----------------------------------------------------------------------
!-- Solves A*X = B for a vector X, where A is specified by the array
!-- u, w, v as returned by svdcmp. m and n are the logical dimensions
!-- of a, and will be equal for square matrices. mp and np are the
!-- physical dimensions of a. b(1:m) is the input right-hand side.
!-- x(1:n) is the output solution vector. No input quantities are
!-- destroyed, so the routine may be called sequentially with 
!-- different b's.
!----------------------------------------------------------------------

      implicit none

      integer, parameter :: dbl=selected_real_kind(12)

      integer m,mp,n,np,nmax
      real(kind=dbl) b(mp),u(mp,np),v(np,np),w(np+1),x(np)
      parameter (nmax = 500)    !Maximum anticipated value of n.
      integer i,j,jj
      real(kind=dbl) s,tmp(nmax)
      do j = 1, n               !Calculate U(T)*B
         s = 0.0
         if (w(j).ne.0.0) then  !Nonzero result only if w(j) is nonzero.
            do i = 1, m
               s = s + u(i,j)*b(i)
            enddo
            s = s/w(j)          !This is the divide by w(j).
         endif
         tmp(j) = s
      enddo
      do j = 1, n               !Matrix multiply by V to get answer.
         s = 0.0
         do jj = 1, n
            s = s + v(j,jj)*tmp(jj)
         enddo
         x(j) = s
      enddo
      return
      end











