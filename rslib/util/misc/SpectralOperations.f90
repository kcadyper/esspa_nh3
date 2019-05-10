!<f90File>**************************************************************
!
! CONTACT:
!
!   Atmospheric & Environmental Research, Inc
!   131 Hartwell Ave
!   Lexington ,MA 02421-3126 USA
!   Phone: 781.761.2288
!   E-mail: guymin@aer.com
!
! COPYRIGHT NOTICE:
!
!   Copyright AER, Inc 2001-2009, All Rights Reserved
!   See the file README-DATARIGHTS.txt included with this release
!   for additional details.
!
!*************************************************************</f90File>

module SpectralOperations

  ! Simple operations for handling spectral data

  implicit none
  private 
  public frqConvert,frqConvertGeneral,spectralConvert
  public NO_EXTRAP_SPECT,EXTRAP_SPECT_LINEAR,EXTRAP_SPECT_CONST

  ! These control behavior when the output frequency grid goes to higher frequencies
  ! than the input frequency grid.  The lower frequency case is not affected.
  integer, parameter :: NO_EXTRAP_SPECT     = 1 ! Do not extrapolate
  integer, parameter :: EXTRAP_SPECT_LINEAR = 2 ! Linear extrapolation
  integer, parameter :: EXTRAP_SPECT_CONST  = 3 ! Constant extrapolation

CONTAINS
  
  subroutine frqconvert(MXCHAN,NFRQMAX,nfrqout,nfrqin,frqin,frqout,polout,cvt, &
       extrapOption)

    ! Creates a matrix that can be used to take emissivities on standard
    ! frequencies and V/H polarization and interpolate to the frequencies
    ! and polarizations of a channel set.
    ! Polarizations P,M,L, and R are assumed to be average of V and H.
    ! The standard-frequency data are assumed to be ordered as V1,H1,V2,H2,...
    ! Example use of cvt:
    ! emisOut(1:nfrqout) = matmul(cvt(1:nfrqout,1:2*nfrqin),emisIn(1:2*nfrqin))

    implicit none

    ! Input Variables
    integer                              :: MXCHAN  ! Channels (1st) dimension of cvt
    integer                              :: NFRQMAX ! Half of 2nd dimension of cvt 
                                                    ! (max # of standard freqencies)
    integer                              :: nfrqout ! Number of channels
    integer                              :: nfrqin  ! Number of standard frequencies
    real,    dimension(NfrqMax)          :: frqin   ! Standard frequencies
    real,    dimension(mxchan)           :: frqout  ! Frequencies of channels
    integer, dimension(mxchan)           :: polout  ! Polarizations of channels
                                                    !  0=V,1=H,2=P,3=M,4=R,5=L
    integer, optional                    :: extrapOption

    ! Output Variable
    real,    dimension(MXCHAN,NFRQMAX*2) :: cvt     ! Conversion matrix

    ! Local Variables
    integer                              :: k, ih, ih1,ih2
    real, dimension(NfrqMax)             :: frqinc
    real                                 :: wi1,wi2,wpV,wpH
    logical                              :: pastInput
    integer                              :: myExtrapOption

    if (present(extrapOption)) then
      myExtrapOption = extrapOption
    else
      myExtrapOption = NO_EXTRAP_SPECT
    endif

    if (myExtrapOption /= NO_EXTRAP_SPECT .and. &
        myExtrapOption /= EXTRAP_SPECT_LINEAR .and. &
        myExtrapOption /= EXTRAP_SPECT_CONST) then
      print *,'err[frqconvert] Invalid extrapOption: ',myExtrapOption
      call errorHalt(1)
    endif

    cvt = 0.

    ! Compute input frequency deltas
    do ih = 2,nfrqin
       frqinc(ih) = frqin(ih) - frqin(ih-1)
    enddo

    ! Loop over all output frequencies
    do k = 1,nfrqout

       ! Find bracketing input frequencies
       if (frqout(k).lt.frqin(1)) then
          print *,'err[frqconvert] frqout(k)<frqin(1);',k,frqout(k),frqin(1)
          call errorHalt(1)
       endif
       
       pastInput = .true.
       do ih = 2,nfrqin
          if (frqout(k).le.frqin(ih)) then
            ih1 = ih-1
            ih2 = ih
            pastInput = .false.
            exit
          endif
       enddo
       if (pastInput) then
         if (myExtrapOption == NO_EXTRAP_SPECT) then
           print *,'err[frqconvert] frqout(k)>frqin(nfrqin);', &
             k,nfrqin,frqout(k),frqin(nfrqin)
           call errorHalt(1)
         else
           ih1=nfrqin-1
           ih2=nfrqin
         endif
       endif
       
       if (pastInput .and. myExtrapOption == EXTRAP_SPECT_CONST) then
         wi1 = 0.
         wi2 = 1.
       else
         wi1 = (frqin(ih2)-frqout(k))/frqinc(ih2)
         wi2 = 1. - wi1
       endif         
       
       if (polout(k).eq.0) then
          wpV = 1.
          wpH = 0.
       elseif(polout(k).eq.1) then
          wpV = 0.
          wpH = 1.
       elseif(polout(k).ge.2 .and. polout(k).le.5) then
          wpV = 0.5
          wpH = 0.5
       else
          print *,'err[frqconvert] Invalid polout(k);',k,polout(k)
          call errorHalt(1)
       endif
       
       cvt(k,ih1*2-1) = wpV*wi1
       cvt(k,ih1*2) = wpH*wi1
       cvt(k,ih2*2-1) = wpV*wi2
       cvt(k,ih2*2) = wpH*wi2
    enddo

    return
  end subroutine frqconvert

  subroutine frqConvertGeneral(frqIn,polIn,frqOut,polOut,cvt)
    
    ! Creates a matrix that can be used to take emissivities on a set of
    ! frequencies and polarizations and interpolate to the frequencies 
    ! and polarizations of a channel set (polarizations P,M,L, and R are 
    ! assumed to be equal to V. This is a generalized version of frqConvert
    ! (no ordering of the V/H data is assumed, although V and H separately
    ! assumed to be ordered monotonically).
    ! Example use of cvt:
    ! emisOut(1:nFrqOut) = matmul(cvt(1:nFrqOut,1:nFrqIn),emisIn(1:nFrqIn))
    
    implicit none
    
    ! I/O Variables
    integer, dimension(:),   intent(in)    :: polIn,polOut ! Polarizations of standard 
                                                           ! frequencies and of channels
    real,    dimension(:),   intent(in)    :: frqIn,frqOut ! Standard frequencies 
                                                           ! and frequencies of channels
                                                           ! 0=V,1=H,2=P,3=M,4=R,5=L
    real,    dimension(:,:), intent(inout) :: cvt          ! Conversion matrix
    
    ! Local Variables
    real, parameter                          :: highFrqExtrp=70.
    integer                                  :: nFrqInV, nFrqInH, k, ih
    integer, dimension(size(frqIn))          :: indxInV,indxInH
    real,    dimension(size(frqIn))          :: frqInV,frqInH
    
    if(size(frqIn) /= size(polIn))then
       print *,'err[frqConvertGeneral] size of frqIn /= size of polIn;', &
            size(frqIn),size(polIn)
       call errorHalt(1)
    endif
    
    if(size(frqOut) /= size(polOut))then
       print *,'err[frqConvertGeneral] size of frqOut /= size of polOut;', &
            size(frqOut),size(polOut)
       call errorHalt(1)
    endif
    
    !--V and H Channels
    nFrqInV = 0
    nFrqInH = 0
    indxInV=0
    indxInH=0
    do ih = 1,size(frqIn)
       if(polIn(ih) == 0)then
          nFrqInV = nFrqInV + 1
          indxInV(nFrqInV) = ih
          frqInV(nFrqInV) = frqin(ih)
       elseif(polIn(ih) == 1)then
          nFrqInH = nFrqInH + 1
          indxInH(nFrqInH) = ih
          frqInH(nFrqInH) = frqin(ih)
       endif
    enddo
    if(nFrqInV <= 1 .or. nFrqInH <= 1)then
       print *,'err[frqConvertGeneral] Not Enough Standard Frequencies for Interpolation;', &
            nFrqInV,nFrqInH
       call errorHalt(1)
    endif
    
    cvt = 0.
    
    !--Loop over output frequencies
    do k = 1,size(frqOut)
       
       if (polOut(k) /= 1)then !--V and polarimetric
          
          ! Find bracketing input frequencies of proper polarization
          if (frqOut(k) < frqInV(1)) then
             print *,'err[frqConvertGeneral] frqOut(k) < frqInV(1);', &
                  k,frqOut(k),frqInV(1)
             call errorHalt(1)
          elseif (frqOut(k) > frqInV(nFrqInV)) then
             print *,'msg[frqConvertGeneral] frqOut(k) > frqInV(nFrqInV);', &
                  k,frqOut(k),frqInV(nFrqInV)
             if(frqInV(nFrqInV) < highFrqExtrp) then
                print *,'err[frqConvertGeneral] frqInV(nFrqInV) < highFrqExtrp;', &
                     frqInV(nFrqInV),highFrqExtrp
                call errorHalt(1)
             else
                print *,'extrapolating'
                ih = nFrqInV
             endif
          else
             do ih = 2,nFrqInV
                if (frqOut(k) <= frqInV(ih))exit
             enddo
          endif
          
          cvt(k,indxInV(ih)) = (frqOut(k)-frqInV(ih-1))/(frqInV(ih)-frqInV(ih-1))
          cvt(k,indxInV(ih-1))   = 1. - cvt(k,indxInV(ih))
          
       else !--H
          
          ! Find bracketing input frequencies of proper polarization
          if (frqOut(k) < frqInH(1)) then
             print *,'err[frqConvertGeneral] frqOut(k) < frqInH(1);', &
                  k,frqOut(k),frqInH(1)
             call errorHalt(1)
          elseif (frqOut(k) > frqInH(nFrqInH)) then
             print *,'err[frqConvertGeneral] frqOut(k) > frqInH(nFrqInH);', &
                  k,frqOut(k),frqInH(nFrqInH)
             call errorHalt(1)
          else
             do ih = 2,nFrqInH
                if (frqOut(k) <= frqInH(ih))exit
             enddo
          endif
          
          cvt(k,indxInH(ih)) = (frqOut(k)-frqInH(ih-1))/(frqInH(ih)-frqInH(ih-1))
          cvt(k,indxInH(ih-1))   = 1. - cvt(k,indxInH(ih))
          
       endif
       
    enddo
    
    return
  end subroutine frqConvertGeneral

  subroutine spectralConvert(spectIn,spectOut,cvt,extrapOK)
    
    ! Creates a matrix that can be used to take data on a set of
    ! spectral locations and interpolate to the locations 
    ! of a channel set.
    ! This does not handle polarization, so should be used only when
    ! polarization is not an issue.
    ! Example use of cvt:
    ! emisOut(1:nFrqOut) = matmul(cvt(1:nFrqOut,1:nFrqIn),emisIn(1:nFrqIn))
    
    implicit none
    
    ! I/O Variables
    real,    dimension(:),   intent(in)    :: spectIn  ! Spectral locations of input
    real,    dimension(:),   intent(in)    :: spectOut ! Spectral locations of output
    real,    dimension(:,:), intent(inout) :: cvt          ! Conversion matrix
    logical,                 intent(in),   optional :: extrapOK  ! no error if 
                                                ! output is out of input range
    
    ! Local Variables
    logical :: myExtrapOK,ascend
    integer, dimension(2)                    :: dimSize
    integer                                  :: nIn,nOut
    integer                                  :: ii, ii1, ii2, k
    
    if (present(extrapOK)) then
      myExtrapOK=extrapOK
    else
      myExtrapOK=.FALSE.
    endif
   
    nIn=size(spectIn)
    nOut=size(spectOut)
    dimSize=shape(cvt)
    if (dimSize(1) < nOut .or. dimSize(2) < nIn) then
       print *,'err[SpectralOperations::spectralConvert] cvt not big enough;', &
            'shape(cvt),nOut,nIn:',dimSize,nOut,nIn
       call errorHalt(1)
    endif
    
    if (all((spectIn(2:nIn)-spectIn(1:nIn-1)) > 0.)) then
      ascend=.TRUE.
    elseif (all((spectIn(2:nIn)-spectIn(1:nIn-1)) < 0.)) then
      ascend=.FALSE.
    else
      print *,'err:[SpectralOperations::spectralConvert]: Spectral locations'
      print *,'not sequential;',spectIn
      call errorHalt(1)
    endif
    
    if (.not. myExtrapOK) then
      if (any(spectOut > maxval(spectIn)) .or. &
          any(spectOut < minval(spectIn))) then
        print *,'err:[SpectralOperations::spectralConvert]: ', &
          'Extrapolation disallowed'
        print *,'maxval(spectOut),maxval(spectIn),', &
          'minval(spectOut),minval(spectIn):', &
          maxval(spectOut),maxval(spectIn),minval(spectOut),minval(spectIn)
        call errorHalt(1)
      endif
    endif

    cvt = 0.

    !--Loop over output spectral locations
    do k = 1,nOut
      if (ascend) then
        do ii=2,nIn
          if (spectIn(ii) >= spectOut(k)) exit
        enddo
      else
        do ii=2,nIn
          if (spectIn(ii) <= spectOut(k)) exit
        enddo
      endif
      ii2=min(ii,nIn)
      ii1=ii2-1
      cvt(k,ii1)=(spectIn(ii2)-spectOut(k))/(spectIn(ii2)-spectIn(ii1))
      cvt(k,ii2)=1.-cvt(k,ii1)
    enddo
    
    return
  end subroutine spectralConvert

end module SpectralOperations
