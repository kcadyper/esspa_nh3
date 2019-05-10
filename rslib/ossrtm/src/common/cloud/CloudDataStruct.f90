MODULE CloudDataStruct
!
! Data structure module used by both cloud scene, and cloud
! optical table for computing cloud optical properties such as
! absorption, scattering, asymmetry, etc.
!
! AER, 01/05/2016
!
  USE CloudParameters
  use params_module, ONLY: mxlay

  IMPLICIT NONE

  PUBLIC

  real, dimension(:), allocatable, private :: spectralPoints

  ! Definition of a physical property at one point in space
  TYPE CloudPhysicalDescr
     character (len=mxNameLength) :: name
     character (len=mxNameLength) :: units
  END TYPE CloudPhysicalDescr
  
  ! Collection of all physical properties (primary size, secondary
  ! size, temperature, melt, density, etc.), at a single point in
  ! space of a hydrometeor
  TYPE CloudPhysicalSet
     integer :: nProperty
     real, allocatable, dimension(:) :: value
  END TYPE CloudPhysicalSet

  ! Definition of a single cloud profile of a hydrometeor
  TYPE CloudScene
     character (len=mxNameLength) :: name
     integer :: nLevel
     integer :: nLayer
     integer :: nProperty
     type (CloudPhysicalDescr), dimension(:), allocatable    :: physDescr
     Type (CloudPhysicalSet), dimension(:), allocatable :: physProp
     Type (CloudPhysicalSet), dimension(:), allocatable :: physPropLayer
     real, dimension(:), allocatable :: hydrAmt
  END TYPE CloudScene

  ! Definition of the grid of a physical property, where the cloud
  ! optical table of a hydrometeor is provided
  TYPE CloudTableGrid
     integer :: id 
     integer :: size
     character (len=mxNameLength) :: name
     character (len=mxNameLength) :: units
     character (len=mxNameLength) :: intplMethod
     integer :: unifGrid
     real :: step
     real, dimension(:), allocatable :: grid
  END TYPE CloudTableGrid

  ! Definition of the cloud optical table of a single cloud optical
  ! property such as absorption, scattering, asymmetry, etc, with the
  ! maximum number of dependencies in terms of the physical properties
  ! plus one (for spectral dimension)
  TYPE CloudOpticalProperty
     integer :: id 
     character (len=mxNameLength) :: name
     character (len=mxNameLength) :: units
     real, dimension(:), allocatable :: values
     real, dimension(:,:), allocatable :: deriv
  END TYPE CloudOpticalProperty

  ! Collection of all optical properties (absorption, scattering,
  ! asymmetry, etc.), at a single point in space of a hydrometeor
  TYPE CloudOpticalSet
     character (len=mxNameLength) :: name
     integer :: nProperty
     TYPE(CloudOpticalProperty), dimension(:), allocatable :: property
  END TYPE CloudOpticalSet

  ! Definition of the cloud optical table of a single cloud optical
  ! property such as absorption, scattering, asymmetry, etc, with the
  ! maximum number of dependencies in terms of the physical properties
  ! plus one (for spectral dimension)
  TYPE CloudOpticalTable
     integer :: id 
     character (len=mxNameLength) :: name
     character (len=mxNameLength) :: units
     real, dimension(:,:), allocatable :: tbl2
     real, dimension(:,:,:), allocatable :: tbl3
     real, dimension(:,:,:,:), allocatable :: tbl4
     real, dimension(:,:,:,:,:), allocatable :: tbl5
     real, dimension(:,:,:,:,:,:), allocatable :: tbl6
  END TYPE CloudOpticalTable

  ! Definition of all optical properties (absorption, scattering,
  ! asymmetry, etc) and the grids of all the independent physical
  ! properties (primary size, secondary size, temperature, density,
  ! melt fraction, spectral points, etc.) for a hydrometeor (liquid,
  ! ice, rain, snow, graupel, melt, etc.)
  TYPE CloudTable
     integer :: id
     integer :: noProperty
     integer :: npProperty
     character (len=mxNameLength) :: name
     TYPE (CloudOpticalTable), dimension(:), allocatable :: oProperty
     TYPE (CloudTableGrid), dimension(:), allocatable :: pProperty
  END TYPE CloudTable

  ! Datatype to handle transition between atmospheric level and
  ! virtual levels inserted due to cloud slabe boundaries

  type, public :: CloudProfileType
     logical                  :: toUse
     integer                  :: levNum
     real, dimension(mxlay+1) :: presLevel
     real, dimension(mxlay+1) :: tempLevel

     integer, dimension(mxlay):: layerIndex
     real, dimension(mxlay)   :: presLayer
     real, dimension(mxlay)   :: tempLayer

     real, dimension(mxlay)   :: molAbs
     real, dimension(mxlay)   :: tau
     real, dimension(mxlay+1) :: ratio

     integer, dimension(2)    :: cloudTop
     integer, dimension(2)    :: cloudBot
     logical, dimension(2)    :: cloudPresent
  end type CloudProfileType

END MODULE CloudDataStruct
