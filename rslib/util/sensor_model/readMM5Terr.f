
      Subroutine readMM5Terr(terrunit,iix,jjx,iixGrid,jjxGrid,
     $     new_num2d,ilat,ilon,new_2dparams)
c     Reads MM5 terrain file. Used by sensor model mm5preproc.
c     Note: Grid read size is 1 element larger than input iixGrid
c     and jjxGrid because we only want coordinates for grid cell center
c     but file gives grid cell corners.
C     2004/11/10:  John Galantowicz
C     Copyright AER, Inc. 2004

      implicit none

      real, dimension(:,:), allocatable :: lat,lon 
      real new_2dparams(iix,jjx,new_num2d)
      character*80 head
      character*10 head2
      integer i,j,cnt,iixGrid,jjxGrid,dum,loop
      integer iix,jjx,new_num2d,ilat,ilon,terrunit

      do i = 1,34 
         read(terrunit,'(A80)')head
      enddo

      allocate(lat(iixGrid+1,jjxGrid+1))
      allocate(lon(iixGrid+1,jjxGrid+1))

      cnt = 0
      loop = 0
      do while (cnt .lt. jjxGrid)    
         do i = 1,2
            read(terrunit,'(A80)')head
         enddo
         do i = iixGrid+1,1,-1
            read(terrunit,*)dum,(lat(i,j+(loop*17)),j= 1,17)
         enddo
         do i = 1,4
            read(terrunit,'(A80)')head
         enddo
         cnt = cnt + 17
         loop = loop + 1
      enddo

      do while (cnt .eq. jjxGrid)    
         do i = 1,2
            read(terrunit,'(A80)')head
         enddo
         do i = iixGrid+1,1,-1
            read(terrunit,*)dum,(lat(i,j+(loop*17)),j= 1,1)
         enddo

         do i = 1,4
            read(terrunit,'(A10)')head2
         enddo
         cnt = cnt + 1
      enddo

c     Now read in lons
      do i = 1,10
         read(terrunit,'(A80)')head
      enddo


      cnt = 0
      loop = 0
      do while (cnt .lt. jjxGrid)    
         do i = 1,2
            read(terrunit,'(A80)')head
         enddo
         do i = iixGrid+1,1,-1
          read(terrunit,'(I3,17(f7.0))')dum,(lon(i,j+(loop*17)),j= 1,17)
         enddo
         do i = 1,4
            read(terrunit,'(A80)')head
         enddo
         cnt = cnt + 17
         loop = loop + 1
      enddo

      do while (cnt .eq. jjxGrid)    
         do i = 1,2
            read(terrunit,'(A80)')head
         enddo
         do i = iixGrid+1,1,-1
            read(terrunit,'(I3,1(f7.0))')dum,(lon(i,j+(loop*17)),j= 1,1)
         enddo

         do i = 1,4
            read(terrunit,'(A10)')head2
         enddo
         cnt = cnt + 1
      enddo

c     Move to output grid (may be smaller) and divide by the 
c     scale factor to get degrees.

      do i = 1,iix
         do j = 1,jjx
            new_2dparams(i,j,ilat) = lat(i,j)/1.0e3
            new_2dparams(i,j,ilon) = lon(i,j)/1.0e2
         enddo
      enddo

      return
      end
