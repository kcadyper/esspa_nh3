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

module gridProf

! <f90Module>***********************************************************
!
! NAME:
!
!   gridProf
!
! PURPOSE:
!
!   Tools for vertical grid remapping of profile data
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  implicit none
  private
  public :: setGridProf,regridProf,closeGridProf

  ! regression matrix and offset
  real,     dimension(:,:), allocatable  :: coef
  real,     dimension(:), allocatable    :: coef0

  ! temp prof vectors
  real,     dimension(:), allocatable    :: profx,profy,prof

  integer :: itop,ibot
  real, dimension(:), allocatable :: pstdLog,xpl

  ! Input and output pressure vectors and number of levels
  type Regres
     real,     dimension(:), pointer    :: prefout,prefin
     integer                            :: nlevout,nlevin
  end type Regres
  type(Regres) RTemp,RH2o,RO3


CONTAINS
  subroutine setGridProf(regrFile,pstd,zone_thk,qCutOffPres)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   setGridProf
!
! PURPOSE:
!
!   Read coefficients for interpolation/extrapolation to output grid
!
! SYNTAX:
!
!   CALL setGridProf(regrFile, pstd, zone_thk, qCutOffPres)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   regrFile     CHAR  path of regression file
!   zone_thk*    REAL  Ln(lower-p/upper-p) of transition zone
!   qCutOffPres* REAL  pressure level at bottom of transition zone
!   
!   INPUTS/OUTPUTS:
!   
!   pstd         REAL  pressure grid
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

!     DESCRIPTION: READ COEFFICIENTS 
!     FOR INTERPOLATION/EXTRAPOLATION TO OUTPUT GRID
!     
!     regrFile    character(len=256) name of regression file
!     pstd        pressure grid
!     zone_thk    ln(lower-p/upper-p) of transition zone
!     qCutOffPres pressure level of bottom of transition zone
!                   default is top of model grid
    use ncdf_module
    implicit none
    character(len=256), intent(in)           :: regrFile
    real, dimension(:), pointer              :: pstd
    real,               intent(in), optional :: zone_thk
    real,               intent(in), optional :: qCutOffPres
    integer :: iurgr
    integer :: i,szP
    real :: h2obot, h2otop
    real :: eps   = 1.0e-2  ! tolerance = 1/100 of pressure
    real :: zone_thk_default

    ! Set default value of transition zone thickness
    zone_thk_default = alog(100./50.)

    ! Open and read regression coefs
    call openNcdfFile(iurgr,RegrFile,'unknown')

    RTemp%nlevout = readNcdfDim(iurgr,'t_nlev_out')
    RTemp%nlevin =  readNcdfDim(iurgr,'t_nlev_in')

    allocate(RTemp%prefout(RTemp%nlevout))
    allocate(RTemp%prefin(RTemp%nlevin))
    call readNcdfData(iurgr,RTemp%prefout,'Tprefout')
    call readNcdfData(iurgr,RTemp%prefin,'Tprefin')

    RH2o%nlevout = readNcdfDim(iurgr,'h2o_nlev_out')
    RH2o%nlevin =  readNcdfDim(iurgr,'h2o_nlev_in')

    allocate(RH2o%prefout(RH2o%nlevout))
    allocate(RH2o%prefin(RH2o%nlevin))
    call readNcdfData(iurgr,RH2o%prefout,'H2oprefout')
    call readNcdfData(iurgr,RH2o%prefin,'H2oprefin')

    RO3%nlevout = readNcdfDim(iurgr,'o3_nlev_out')
    RO3%nlevin =  readNcdfDim(iurgr,'o3_nlev_in')
    if (RO3%nlevout >0 .and. RO3%nlevin>0) then
       allocate(RO3%prefout(RO3%nlevout))
       allocate(RO3%prefin(RO3%nlevin))
       call readNcdfData(iurgr,RO3%prefout,'O3prefout')
       call readNcdfData(iurgr,RO3%prefin,'O3prefin')
    endif

    allocate(coef0(RH2o%nlevout+RTemp%nlevout+RO3%nlevout))
    allocate(coef(RTemp%nlevout+RH2o%nlevout+RO3%nlevout,&
         RTemp%nlevin+RH2o%nlevin+RO3%nlevin))
    call readNcdfData(iurgr,coef0,'coef0')
    call readNcdfData(iurgr,coef,'coef')

    call closeNcdfFile(iurgr)

    ! Set size of temp arrays
    allocate(prof(RH2o%nlevout+RTemp%nlevout+RO3%nlevout))
    allocate(profx(RH2o%nlevin+RTemp%nlevin+RO3%nlevin))
    allocate(profy(RH2o%nlevin+RTemp%nlevin+RO3%nlevin))

    ! Do size check
    szP = size(RTemp%prefout)

    ! Top if input grid determines bottom of transition zone for cubic interpolation
    
    allocate(pstd(size(RH2o%prefout)),pstdLog(size(RH2o%prefout)))
    allocate(xpl(size(RH2o%prefout)))

    pstd = RH2o%prefout
    pstdLog = alog(pstd)
    h2obot = alog(RH2o%prefin(1))
    if (present(qCutOffPres)) then
       if (qCutOffPres > 0.) h2obot = alog(qCutOffPres)  ! negative => default
    endif
    
    do i=1,szP ! Determine indices (in terms of output grid) of pressure levels 
               !   just above input grid
       if(pstdLog(i).ge.(h2obot-eps))then
          ibot = i 
          exit
       endif
       if(i.eq.szP)then
          print *,'err[gridProf::setGridProf]: ', &
            'Could not determine grid index for bottom of transition zone'
          call errorHalt(1)
       endif
    enddo

    ! If transition zone not specified, default to zone of zero depth
    if(present(zone_thk))then
       h2otop = pstdLog(ibot) - zone_thk
       print*,'Using input zone_thk=',zone_thk
    else
       h2otop = pstdLog(ibot) - zone_thk_default
       print*,'Using default zone_thk=',zone_thk_default
    endif

    do i=1,szP ! Determine indices of pressure levels just above a transition zone
       if(pstdLog(i).gt.h2otop)then
          itop = i
          exit
       endif
    enddo

    print*,'Bottom and top of transition zone are (mb)', &
       ibot,itop,pstd(ibot),pstd(itop)

  end subroutine setGridProf

  subroutine regridProf(profGrid,IGgrid,NGgrid,MolIDgrid, &
     profFov,IGfov,NGfov,T2m,q2m)


!<f90Subroutine>********************************************************
!
! NAME:
!
!   regridProf
!
! PURPOSE:
!
!   Execute interpolation/extrapolation of a profile to output grid
!
! SYNTAX:
!
!   CALL regridProf(profGrid, IGgrid, NGgrid, MolIDgrid, profFov, IGfov, NGfov, 
!      T2m, q2m)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   profGrid  REAL                Atmospheric profile on grid
!   IGgrid    TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                 geophysical state vector on grid 
!   NGgrid    TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                 geophysical state vector on grid 
!   IGfov     TYPE(STATEINDEX_T)  Starting indices for sections of 
!                                 geophysical state vector at FOV 
!   NGfov     TYPE(STATEINDEX_T)  Number of elements for sections of 
!                                 geophysical state vector at FOV 
!   MolIDgrid INTEGER             List of IDs of relevant molecular species 
!                                 at FOV 
!   T2m*      REAL                2-m temperature
!   q2m*      REAL                2-m humidity
!   
!   INPUTS/OUTPUTS:
!   
!   profFov   REAL                Atmospheric profile at FOV
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    use StateIndexModule
    use constants, only: g0,Rd
    implicit none
    real, dimension(:),   intent(in)            :: profGrid
    type(StateIndex_t),   intent(in)            :: IGgrid, NGgrid
    real, dimension(:),   intent(inout)         :: profFov
    type(StateIndex_t),   intent(in)            :: IGfov, NGfov
    integer, dimension(:),intent(in)            :: MolIDgrid
    real,                 intent(in),  optional :: T2m,q2m
    real :: q2mint, T2mint, q2mext, T2mext, q2madj, T2madj
    real :: Tavg, pavg, rodens, dp, dz, dTdz, qavg, dqdz
    integer :: i,j,k,ilev
    real :: a,svp
    real :: dx,y1,y2,g1,g2,g3,b,c,d,y
    real :: c_fit,l_fit,l_lo,l_hi
    real, parameter :: factor=0.05     ! fraction of ln(qv) used to determine tolerance of cubic fit
    real, parameter :: h2o_min=1.e-12  ! smallest number to represent moisture (due alog function)
    real, parameter :: max_dtdz=0.0050 ! abs(largest permiss T gradient in downward extrap) 
                                       ! currently set to 3x dry adiabatic lapse rate
    real, parameter :: max_dqdz=1.e-8  ! abs(largest permiss q gradient in downward extrap)
                                       ! currently set to approx. saturation MR lapse rate
    real, parameter :: onePlusFactor=1.+factor,oneMinusFactor=1.-factor
    real :: ln_h2o_min                 ! ln of h2o_min
    integer :: nMol, h2oid,O3ID
    integer :: H2oIDgrid,O3IDgrid,H2oIDfov,O3IDfov
    
! Figure out the H2O ID
    nMol = getNMol(NGfov)
    H2oID = whereH2O(MolIDgrid)
    O3ID = searchID(molIDgrid,idO3)

    !if (nc1 /= NGgrid%Temp+NGgrid%mol(H2Oid)) then
    if (RTemp%nlevin+RH2o%nlevin /= NGgrid%Temp+NGgrid%mol(H2Oid)) then
         print *,"err[gridProf::regridProf]:  inconsistent #mapping grids."
       call errorHalt(1)
    end if
    if(O3ID>0) then
       if (RO3%nlevin /= NGgrid%mol(O3ID)) then
          print *,"err[gridProf::regridProf]:  inconsistent #mapping O3 grids."
          call errorHalt(1)
       endif
    endif

    ln_h2o_min = alog(h2o_min)

    ! Map input prof to profx, profy vector; determine ordering of profile
    profx=0.
    profx(1:NGgrid%Temp) = profGrid(IGgrid%Temp:IGgrid%Temp+NGgrid%Temp-1)
    profx(IGgrid%Temp+NGgrid%Temp:&
         IGgrid%Temp+NGgrid%Temp+NGgrid%mol(H2Oid)-1) = &
         max(profGrid(IGgrid%mol(H2Oid):IGgrid%mol(H2Oid)+NGgrid%mol(H2Oid)-1),ln_h2o_min)
    if(O3ID>0) profx(IGgrid%Temp+NGgrid%Temp+NGgrid%mol(H2Oid):&
         IGgrid%Temp+NGgrid%Temp+NGgrid%mol(H2Oid)+NGgrid%mol(O3ID)-1) = &
         max(profGrid(IGgrid%mol(O3ID):IGgrid%mol(O3ID)+NGgrid%mol(O3ID)-1),ln_h2o_min)

    profy=0.
    profy(1:NGgrid%Temp) = profGrid(IGgrid%Temp:IGgrid%Temp+NGgrid%Temp-1)
    profy(IGgrid%Temp+NGgrid%Temp:&
         IGgrid%Temp+NGgrid%Temp+NGgrid%cldliq-1) = &
         alog(max(profGrid(IGgrid%cldliq:IGgrid%cldliq+NGgrid%cldliq-1),h2o_min))

! Interpolate/extrapolate to other levels
    ! Map input grid to output grid(only if good prof)
    if(profx(1).eq.-9999) then
       prof(1:NGfov%Temp+NGfov%mol(H2Oid))=-9999
    else
       prof=matmul(coef,profx)
       prof=prof+coef0
! Cubic interpolation of moisture for points within transition zone computed 
! in subroutine SetGridProf. Algorithm from A. Lipton (June 2005)
       dx = pstdLog(ibot)-pstdLog(itop)
       xpl(1:ibot-itop+1) = pstdLog(itop:ibot)-pstdLog(itop)
       y1 = prof(NGfov%Temp+itop)
       y2 = prof(NGfov%Temp+ibot)
       g1 = (prof(NGfov%Temp+itop)-prof(NGfov%Temp+itop-1))/ &
         (pstdLog(itop)-pstdLog(itop-1))
       g2 = (prof(NGfov%Temp+ibot+1)-prof(NGfov%Temp+ibot))/ &
         (pstdLog(ibot+1)-pstdLog(ibot))
       g3 = (prof(NGfov%Temp+ibot)-prof(NGfov%Temp+itop))/dx
       a = y1
       b = g1
       c = (3*(y2-y1)-(2*g1+g2)*dx)/dx**2
       d = ((g1+g2)*dx+2*(y1-y2))/dx**3
       ! Compute cubic fit and error bounds derived from a linear fit between ibot and itop
       do i=1,ibot-itop+1     !-1
          c_fit = a+b*xpl(i)+c*xpl(i)**2+d*xpl(i)**3
          l_fit = a+g3*xpl(i)
          l_lo = min(l_fit*onePlusFactor,l_fit*oneMinusFactor)   ! l_fit may be + or -
          l_hi = max(l_fit*onePlusFactor,l_fit*oneMinusFactor)   ! l_fit may be + or -
          prof(NGfov%Temp+itop-1+i) = min(l_hi,max(l_lo,c_fit))
       enddo
       prof(NGfov%Temp+1:NGfov%Temp+NGfov%mol(H2Oid))= &
            exp(prof(NGfov%Temp+1:NGfov%Temp+NGfov%mol(H2Oid)))
    endif

! Fill in output vector. It is assumed that the NWP vector consists of 
! T and H2O, followed by Tskin, Psfc, CLW, and wind.

    profFov(IGfov%Temp:IGfov%Temp+NGfov%Temp-1)=prof(1:NGfov%Temp)
    profFov(IGfov%Tskin)=profGrid(IGgrid%Tskin)
    profFov(IGfov%Psfc)=profGrid(IGgrid%Psfc)
    profFov(IGfov%mol(H2Oid):IGfov%mol(H2Oid)+NGfov%mol(H2Oid)-1)= &
         prof(NGfov%Temp+1:NGfov%Temp+NGfov%mol(H2Oid))
    if(O3ID>0) profFov(IGfov%mol(O3ID):IGfov%mol(O3ID)+NGfov%mol(O3ID)-1)= &
         exp(prof(NGfov%Temp+NGfov%mol(H2Oid)+1:NGfov%Temp+NGfov%mol(H2Oid)+ &
         NGfov%mol(O3ID)))
    profFov(IGfov%wind:IGfov%wind+NGfov%wind-1) = &
         profGrid(IGgrid%wind:IGgrid%wind+NGgrid%wind-1)    

   !  Repeat for CLW, but do not extrapolate 
    if(NGgrid%cldliq .gt. 0) then
       if(profy(1).eq.-9999) then
          prof(1:NGfov%Temp+NGfov%mol(H2Oid))=-9999
       else
          prof=matmul(coef,profy)
          prof=prof+coef0
          prof(NGfov%Temp+1:NGfov%Temp+NGfov%cldliq)= &
               exp(prof(NGfov%Temp+1:NGfov%Temp+NGfov%cldliq))
          !=================================================
          !--Determine index of first pressure level above surface:
          !   This assumes CLW profile has same levels as H2O profile.
          do ilev=size(RH2o%prefout),1,-1 
             if (RH2o%prefout(ilev).lt.profGrid(IGgrid%Psfc)) exit
          enddo
          !--Do not allow subsurface CLW extrapolation to exceed first
          !    CLW value above the surface:
          do i=(NGfov%Temp+NGfov%cldliq),(NGfov%Temp+(ilev)),-1
             prof(i) = min(prof(i),prof(NGfov%temp+ilev))
          enddo
          !=================================================
       endif
    endif

! If CLW is present, fill in output vector; zero-out above input profile (i.e., ibot)
    if(NGgrid%cldliq .gt. 0) then
       profFov(IGfov%cldliq:IGfov%cldliq+NGfov%cldliq-1)= &
            max(h2o_min,prof(NGfov%Temp+1:NGfov%Temp+NGfov%cldliq))
       profFov(IGfov%cldliq:IGfov%cldliq-1+ibot-1) = h2o_min
    endif

! If 2-m temperature or humidity available, use to ensure 
    if(present(T2m))then
       ! Find output grid index of first below ground p-level
       if(profFov(IGfov%psfc) .le. RTemp%prefout(RTemp%nlevout)) then
          i=1
          do while(RTemp%prefout(i) .lt. profFov(IGfov%psfc) .and. i .le. RTemp%nlevout-1)
             i=i+1
          enddo
          !Temperature
          !Compute temperature lapse rate, dT/dz; limit to 
          !pre-specified value
          Tavg = (T2m + profFov(IGfov%Temp+(i-1)-1))*0.5
          pavg = (profFov(IGfov%psfc) + RTemp%prefout(i-1))*0.5
          rodens = 100.*pavg/(Rd*Tavg)
          dp = (profFov(IGfov%psfc)-RTemp%prefout(i-1))*100.
          dz = dp/(-rodens*g0)
          dTdz = (profFov(IGfov%Temp+(i-1)-1) - T2m)/dz
          T2madj = profFov(IGfov%Temp+(i-1)-1) - sign(1.,dTdz)*min(abs(dTdz),max_dTdz)*dz
          do k=i,RTemp%nlevout
             call lint(RTemp%prefout(i-1),profFov(IGfov%psfc), &
                  profFov(IGfov%Temp+(i-1)-1),T2madj, &
                  RTemp%prefout(k),T2mext)
             profFov(IGfov%Temp+(k-1)) = T2mext
             call lint(RTemp%prefout(i-1),RTemp%prefout(k), &
                  profFov(IGfov%Temp+(i-1)-1),T2mext, &
                  profFov(IGfov%psfc),T2mint)
             if(T2mint-T2madj .gt. 0.1) print *, &
                  'WARNING: T2mint .ne. T2madj',T2mint,T2madj
          enddo
       endif
    endif
    if(present(q2m))then
       if(profFov(IGfov%psfc) .le. RH2o%prefout(RH2o%nlevout)) then
          i=1
          do while(RH2o%prefout(i) .lt. profFov(IGfov%psfc) .and. i .le. RH2o%nlevout-1)
             i=i+1
          enddo
          !Specific humidity
          !Compute sh lapse rate, dqv/dz; limit to
          !pre-specified value
          qavg = (q2m + profFov(IGfov%mol(H2Oid)+(i-1)-1))*0.5
          pavg = (profFov(IGfov%psfc) + RH2o%prefout(i-1))*0.5
          rodens = 100.*pavg/(Rd*Tavg)
          dp = (profFov(IGfov%psfc)-RH2o%prefout(i-1))*100.
          dz = dp/(-rodens*g0)
          dqdz = (profFov(IGfov%mol(H2Oid)+(i-1)-1) - q2m)/dz
          q2madj = profFov(IGfov%mol(H2Oid)+(i-1)-1) - sign(1.,dqdz)*min(abs(dqdz),max_dqdz)*dz
          do k=i,RH2o%nlevout
             call lint(RH2o%prefout(i-1),profFov(IGfov%psfc), &
                  profFov(IGfov%mol(H2Oid)+(i-1)-1),q2madj, &
                  RH2o%prefout(k),q2mext)
             profFov(IGfov%mol(H2Oid)+(k-1)) = q2mext
             call lint(RH2o%prefout(i-1),RH2o%prefout(k), &
                  profFov(IGfov%mol(H2Oid)+(i-1)-1),q2mext, &
                  profFov(IGfov%psfc),q2mint)
             if(q2mint-q2madj .gt. 0.1) print*,'WARNING: q2mint .ne. q2madj'            
          enddo
       endif
    endif

! Optional: check humidities for super-saturation.
! The indices need to be fixed when this is turned on.
!
!      do 230 i=41,61
!         a=svp(profFov(i-21))
!         a=alog(622.*a/(pref(i-21)-a))
!         profFov(i)=min(profFov(i),a)
! 230  continue

    return
  end subroutine regridProf

  subroutine closeGridProf()

!<f90Subroutine>********************************************************
!
! NAME:
!
!   closeGridProf
!
! PURPOSE:
!
!   Clean up.
!
! SYNTAX:
!
!   CALL closeGridProf()
!
! ARGUMENTS:
!
!   
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    implicit none
    deallocate(RTemp%prefout,RTemp%prefin,RH2o%prefout,RH2o%prefin)
    if(associated(RO3%prefout)) deallocate(RO3%prefout)
    if(associated(RO3%prefin))  deallocate(RO3%prefin)
    deallocate(coef0,coef)
    deallocate(prof,profx,profy)
    deallocate(pstdLog,xpl)
    return
  end subroutine closeGridProf

end module gridProf
