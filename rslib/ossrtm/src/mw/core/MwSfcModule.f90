MODULE MwSfcModule

  use EmisTab
  implicit none
  private
  public :: MwSfcInit,computeRefl,MwSfcClose
  public :: nGauss,gaussAng

  integer                                    :: nGauss
  real,    dimension(:),        pointer      :: gaussAng
  real                                       :: wcasp
  real                                       :: PlandMinLd

CONTAINS
  
  SUBROUTINE MwSfcInit(fname,outBndsExit,PlandMinLd_in)
    
    !  Initialize access to tabulated emissivity data.
    
    character(len=*), dimension (:), intent(in)  :: fname
    logical, intent(in)                          :: outBndsExit
    real, intent(in)                             :: PlandMinLd_in
    integer                                      :: indx
    
    do indx=1,size(fname)
       if(indx == 1)then
          call emisTabInit(indx,fname(indx),'emisSpec')
       elseif(indx == 2)then
          call emisTabInit(indx,fname(indx),'emWindComp')
       else
          call emisTabInit(indx,fname(indx),'Reflectivity', &
               nGaussOut=nGauss,gaussAngOut=gaussAng,outBndsExit_in=outBndsExit)
          WCASP=16./float(nGauss)
       endif
    enddo

    PlandMinLd=PlandMinLd_in

    return
    
  END SUBROUTINE MwSfcInit
  
  !--------------------------------------------------------------------------
  
  subroutine computeRefl(wspd,obsang,pland,wrefl,drodwind)
    !---Input variables
    REAL,                  INTENT(IN)  :: wspd,pland
    REAL,    DIMENSION(:), INTENT(IN)  :: obsang
    !---Output variables
    REAL,    DIMENSION(:,:), INTENT(INOUT) :: wrefl,drodwind
    !---Local variables
    INTEGER                                 :: ich,i
    REAL                                    :: wspd0
!!$    REAL,dimension(SIZE(obsang),NGAUSS)     :: rcoeff1  !--Needed for quadratic casp

    wrefl=0.
    drodwind=0.

    if(pland > PlandMinLd)return !--Specular for land

!!$    wcasp=4. !---Uncomment for quadratic casp
    wspd0=max(wspd,wcasp)
    
    do ich=1,SIZE(obsang)
       call getEmisTab(3,intrpScheme=2, &
            listChans=(/(i,i=nGauss*(ich-1)+1,nGauss*(ich-1)+nGauss)/), &
            speed=wspd0,angle=abs(obsAng(ich)), &
            emis=wrefl(ich,1:nGauss), &
            demdspeed=drodwind(ich,1:nGauss))
    enddo
    
!!$    rcoeff1=wrefl !---Uncomment for quadratic casp

    if(wspd < wcasp)then
       !---Uncomment for quadratic casp
!!$       wrefl=wrefl+(wspd-wcasp)*drodwind+ &
!!$            (wspd-wcasp)*(wspd-wcasp)*(wcasp*drodwind-wrefl)/wcasp/wcasp
!!$       drodwind=drodwind+(wspd-wcasp)*(wcasp*drodwind-rcoeff1)/wcasp/wcasp*2.
       !---Comment out for quadratic casp
       drodwind=wrefl/wcasp
       wrefl=wspd/wcasp*wrefl
    endif

    return
  end subroutine computeRefl

  subroutine MwSfcClose()

    call emisTabDestroy(0)
    return

  end subroutine MwSfcClose

END MODULE MwSfcModule
