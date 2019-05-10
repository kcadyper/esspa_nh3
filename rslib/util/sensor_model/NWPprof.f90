module NWPprof
  implicit none
  public  :: nwpProfGet
  interface nwpProfGet
  module procedure getNWPprof
  end interface nwpProfGet
  private
  contains
  include "getNWPprof.f90"
end module NWPprof