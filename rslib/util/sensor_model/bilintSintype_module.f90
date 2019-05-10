module bilintSintype_module

! create defined type needed for bilinear interpolation option

  type bilintSin


! two sets of arrays needed due to: 
! (1) shift of grid values w/r to array locations (i.e. value stored at array(1,1) is valid at grid(0.5,0.5))
! (2) possibility of wraparound at dateline. Required array locations will not be consistent with grid values
! arrays for storing *grid* locations (non-integer values, needed for finite differences)
    real, dimension(:), pointer :: i0floorG, i0ceilG, jfloorG, jceilG
    real, dimension(:), pointer :: i_jfloorG, i_jceilG
! arrays for storing *array* index locations within local tiles (from tile files)
    integer, dimension(:), pointer :: i0floorALocal, i1floorALocal
    integer, dimension(:), pointer :: i0ceilALocal, i1ceilALocal
    integer, dimension(:), pointer :: jfloorALocal, jceilALocal
! arrays for storing cols and rows of corresponding local tiles (to identify tile files)
    integer, dimension(:), pointer :: i0floor_tilecol, i1floor_tilecol
    integer, dimension(:), pointer :: i0ceil_tilecol, i1ceil_tilecol
    integer, dimension(:), pointer :: jfloor_tilerow, jceil_tilerow


 end type bilintSin


end module bilintSintype_module
