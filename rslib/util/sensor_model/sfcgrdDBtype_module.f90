module sfcgrdDBtype_module

! create defined type to contain essential database parameters

  type sfcgrdDB

     integer :: ntilefiles

     character(len=256) :: listfile_tilefiles
     character(len=256), allocatable, dimension(:) :: tilefile

     integer, allocatable, dimension(:) :: tile_colglobal
     integer, allocatable, dimension(:) :: tile_rowglobal
     integer, allocatable, dimension(:) :: tile_ncolglobal
     integer, allocatable, dimension(:) :: tile_nrowglobal
     integer, allocatable, dimension(:) :: tile_ncol
     integer, allocatable, dimension(:) :: tile_nrow
     real, allocatable, dimension(:) :: tile_colGridOffset
     real, allocatable, dimension(:) :: tile_rowGridOffset
     real, allocatable, dimension(:) :: tile_Re, tile_scale

  end type sfcgrdDB


end module sfcgrdDBtype_module
