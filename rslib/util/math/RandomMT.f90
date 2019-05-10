! A C-program for MT19937, with initialization improved 2002/1/26.
! Coded by Takuji Nishimura and Makoto Matsumoto.

! Code converted to Fortran 95 by Josi Rui Faustino de Sousa
! Date: 2002-02-01

! Before using, initialize the state by using init_genrand(seed)
! or init_by_array(init_key).

! This library is free software.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

! Copyright (C) 1997, 2002 Makoto Matsumoto and Takuji Nishimura.
! Any feedback is very welcome.
! http://www.math.keio.ac.jp/matumoto/emt.html
! email: matumoto@math.keio.ac.jp

! Updated March 2006 by Katherine Quinn, AER.  Added generic interfaces, 
! defined a structure that includes the state vector (mt(n)) and 
! state counter (mti), and added subroutines to get and set the state.

module RandomMT

  implicit none

  intrinsic :: bit_size

  private
!!!  public  :: init_genrand, init_by_array
!!!  public  :: genrand_int32, genrand_int31
!!!  public  :: genrand_real1, genrand_real2, genrand_real3, genrand_res53
  public  :: initRandMT, randMT, getStateMT, setStateMT, nextRandMT
  public  :: randMTd, nextRandMTd

  interface initRandMT
    module procedure init_genrand_int
    module procedure init_genrand_intarray
  end interface

  integer,  parameter  :: intg = selected_int_kind( 9 )
  integer,  parameter  :: dobl = selected_real_kind( 15, 307 )

  integer,  parameter  :: wi = intg
  integer,  parameter  :: wr = dobl

  integer,  public,  parameter  :: MT_INT_KIND = intg
  integer,  public,  parameter  :: MT_REAL_KIND = dobl


  ! Period parameters
  integer( kind = wi ), parameter :: n = 624_wi
  integer( kind = wi ), parameter :: m = 397_wi
  integer( kind = wi ), parameter :: hbs = bit_size( n ) / 2_wi
  integer( kind = wi ), parameter :: qbs = hbs / 2_wi
  integer( kind = wi ), parameter :: tbs = 3_wi * qbs

  integer( kind = wi )  :: mt(n)                ! the array for the state vector
  logical( kind = wi )  :: mtinit = .false._wi  ! means mt[N] is not initialized
  integer( kind = wi )  :: mti = n + 1_wi       ! mti==N+1 means mt[N] is not initialized

  type, public  :: StateMT_t
    integer( kind = wi )  :: vec(n) 
    integer( kind = wi )  :: cnt = n + 1_wi
    logical( kind = wi )  :: init = .false._wi
  end type StateMT_t

  contains

  elemental function uiadd( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 + b1
    s2 = a2 + b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )

  end function uiadd

  elemental function uisub( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 - b1
    s2 = a2 - b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )

  end function uisub

  elemental function uimlt( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a0, a1, a2, a3
    integer( kind = wi )  :: b0, b1, b2, b3
    integer( kind = wi )  :: p0, p1, p2, p3

    a0 = ibits( a, 0, qbs )
    a1 = ibits( a, qbs, qbs )
    a2 = ibits( a, hbs, qbs )
    a3 = ibits( a, tbs, qbs )
    b0 = ibits( b, 0, qbs )
    b1 = ibits( b, qbs, qbs )
    b2 = ibits( b, hbs, qbs )
    b3 = ibits( b, tbs, qbs )
    p0 = a0 * b0
    p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
    p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs )
    p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
    c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) )
    c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
    c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) )

  end function uimlt

  ! initializes mt[N] with a seed
  subroutine init_genrand( s )

    implicit none

    intrinsic :: iand, ishft, ieor, ibits

    integer( kind = wi ), intent( in )  :: s

    integer( kind = wi )  :: i, mult_a

    data mult_a /z'6C078965'/

    mtinit = .true._wi
    mt(1) = ibits( s, 0, 32 )
    do i = 2, n, 1
      mt(i) = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
      mt(i) = uimlt( mt(i), mult_a )
      mt(i) = uiadd( mt(i), uisub( i, 1_wi ) )
      ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
      ! In the previous versions, MSBs of the seed affect
      ! only MSBs of the array mt[].
      ! 2002/01/09 modified by Makoto Matsumoto
      mt(i) = ibits( mt(i), 0, 32 )
      ! for >32 bit machines
    end do

  end subroutine init_genrand

  ! initialize by an array with array-length
  ! init_key is the array for initializing keys
  ! key_length is its length
  subroutine init_by_array( init_key )

    implicit none

    intrinsic :: iand, ishft, ieor

    integer( kind = wi ), intent( in )  :: init_key(:)

    integer( kind = wi )  :: i, j, k, tp, key_length
    integer( kind = wi )  :: seed_d, mult_a, mult_b, msb1_d

    data seed_d /z'12BD6AA'/
    data mult_a /z'19660D'/
    data mult_b /z'5D588B65'/
    data msb1_d /z'80000000'/

    key_length = size( init_key, dim = 1 )
    call init_genrand( seed_d )
    i = 2_wi
    j = 1_wi
    do k = max( n, key_length ), 1, -1
      tp = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
      tp = uimlt( tp, mult_a )
      mt(i) = ieor( mt(i), tp )
      mt(i) = uiadd( mt(i), uiadd( init_key(j), uisub( j, 1_wi ) ) ) ! non linear
      mt(i) = ibits( mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
      i = i + 1_wi
      j = j + 1_wi
      if ( i > n ) then
        mt(1) = mt(n)
        i = 2_wi
      end if
      if ( j > key_length) j = 1_wi
    end do
    do k = n-1, 1, -1
      tp = ieor( mt(i-1), ishft( mt(i-1), -30 ) )
      tp = uimlt( tp, mult_b )
      mt(i) = ieor( mt(i), tp )
      mt(i) = uisub( mt(i), uisub( i, 1_wi ) ) ! non linear
      mt(i) = ibits( mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
      i = i + 1_wi
      if ( i > n ) then
        mt(1) = mt(n)
        i = 2_wi
      end if
    end do
    mt(1) = msb1_d ! MSB is 1; assuring non-zero initial array
  end subroutine init_by_array

  ! generates a random number on [0,0xffffffff]-interval
  function genrand_int32( ) result( y )

    implicit none

    intrinsic :: iand, ishft, ior, ieor, btest, ibset, mvbits

    integer( kind = wi )  :: y

    integer( kind = wi )  :: kk
    integer( kind = wi )  :: seed_d, matrix_a, matrix_b, temper_a, temper_b

    data seed_d   /z'5489'/
    data matrix_a /z'9908B0DF'/
    data matrix_b /z'0'/
    data temper_a /z'9D2C5680'/
    data temper_b /z'EFC60000'/

    if ( mti > n ) then ! generate N words at one time
      if ( .not. mtinit ) then ! if init_genrand() has not been called, a default initial seed is used
         call init_genrand( seed_d )
         print *,'Warning[RandomMT::genrand_int32]: state has not been initialized.'
         print *,'Using default seed with init_genrand(iseed)'
      end if
      do kk = 1, n-m, 1
        y = ibits( mt(kk+1), 0, 31 )
        call mvbits( mt(kk), 31, 1, y, 31 )
        if ( btest( y, 0 ) ) then
          mt(kk) = ieor( ieor( mt(kk+m), ishft( y, -1 ) ), matrix_a )
        else
          mt(kk) = ieor( ieor( mt(kk+m), ishft( y, -1 ) ), matrix_b )
        end if
      end do
      do kk = n-m+1, n-1, 1
        y = ibits( mt(kk+1), 0, 31 )
        call mvbits( mt(kk), 31, 1, y, 31 )
        if ( btest( y, 0 ) ) then
          mt(kk) = ieor( ieor( mt(kk+m-n), ishft( y, -1 ) ), matrix_a )
        else
          mt(kk) = ieor( ieor( mt(kk+m-n), ishft( y, -1 ) ), matrix_b )
        end if
      end do
      y = ibits( mt(1), 0, 31 )
      call mvbits( mt(n), 31, 1, y, 31 )
      if ( btest( y, 0 ) ) then
        mt(kk) = ieor( ieor( mt(m), ishft( y, -1 ) ), matrix_a )
      else
        mt(kk) = ieor( ieor( mt(m), ishft( y, -1 ) ), matrix_b )
      end if
      mti = 1_wi
    end if
    y = mt(mti)
    mti = mti + 1_wi
    ! Tempering
    y = ieor( y, ishft( y, -11) )
    y = ieor( y, iand( ishft( y, 7 ), temper_a ) )
    y = ieor( y, iand( ishft( y, 15 ), temper_b ) )
    y = ieor( y, ishft( y, -18 ) )

  end function genrand_int32

  ! generates a random number on [0,0x7fffffff]-interval
  function genrand_int31( ) result( i )

    implicit none

    intrinsic :: ishft

    integer( kind = wi )  :: i

    i = ishft( genrand_int32( ), -1 )

  end function genrand_int31

  ! generates a random number on [0,1]-real-interval
  function genrand_real1( ) result( r )

    implicit none

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = genrand_int32( )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = real( a0, kind = wr ) / 4294967295.0_wr
    r = real( a1, kind = wr ) * ( 65536.0_wr / 4294967295.0_wr ) + r
    ! divided by 2^32-1

  end function genrand_real1

  ! generates a random number on [0,1)-real-interval
  function genrand_real2( ) result( r )

    implicit none

    intrinsic :: ibits

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = genrand_int32( )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = real( a0, kind = wr ) / 4294967296.0_wr
    r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32

  end function genrand_real2

  ! generates a random number on (0,1)-real-interval
  function genrand_real3( ) result( r )

    implicit none

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = genrand_int32( )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = ( real( a0, kind = wr ) + 0.5_wr ) / 4294967296.0_wr
    r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32

  end function genrand_real3

  ! generates a random number on [0,1) with 53-bit resolution
  function genrand_res53( )  result( r )

    implicit none

    intrinsic :: ishft

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a0, a1
    integer( kind = wi )  :: b, b0, b1

    a = ishft( genrand_int32( ), -5 )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    b = ishft( genrand_int32( ), -6 )
    b0 = ibits( b, 0, hbs )
    b1 = ibits( b, hbs, hbs )
    r = real( a1, kind = wr ) / 2048.0_wr
    r = real( a0, kind = wr ) / 134217728.0_wr + r
    r = real( b1, kind = wr ) / 137438953472.0_wr + r
    r = real( b0, kind = wr ) / 9007199254740992.0_wr + r

  end function genrand_res53
  ! These real versions are due to Isaku Wada, 2002/01/09 added

  ! Interfaces to genrand_real2 to give single prec real answer
  function randMT() result(snglr)

    implicit none
    real               :: snglr
    real( kind = wr )  :: r

    r = genrand_real2( )
    snglr = real(r)

  end function randMT

  ! Interfaces to genrand_real2 to give double prec real answer
  function randMTd() result(dblr)

    implicit none
    double precision   :: dblr
    real( kind = wr )  :: r

    r = genrand_real2( )
    dblr = dble(r)

  end function randMTd

  ! Interface to init_genrand subroutine with generic integer input
  subroutine init_genrand_int( intgs, state )

    implicit none
    integer, intent( in ) :: intgs
    type(StateMT_t), intent( inout ), optional :: state
    integer( kind = wi )  :: s

    if (range(intgs) .gt. range(s)) then
      print *,'err[RandomMT::init_genrand_int]: ',&
        'Seed number of digits too big'
      STOP
    end if

    mti = n + 1_wi   !--- reset counter
    s = int(intgs, kind = wi)
    call init_genrand( s )

    if (present(state)) call getStateMT(state)

  end subroutine init_genrand_int

  ! Interface to init_by_array subroutine with generic integer input
  subroutine init_genrand_intarray( intgkey, state )

    implicit none
    integer, dimension(:), intent( in )             :: intgkey
    type(StateMT_t), intent( inout ), optional      :: state
    integer( kind = wi ), dimension(size(intgkey))  :: init_key

    if (range(intgkey(1)) .gt. range(init_key(1))) then
      print *,'err[RandomMT::init_genrand_intarray]: ',&
        'Seed number of digits too big'
      STOP
    end if

    mti = n + 1_wi   !--- reset counter
    init_key = int(intgkey, kind = wi)
    call init_by_array( init_key )

    if (present(state)) call getStateMT(state)

  end subroutine init_genrand_intarray

  ! Get the current state vector and counter
  subroutine getStateMT( State0 )

    implicit none
    type(StateMT_t), intent( inout ) :: State0

    State0%vec = mt
    State0%cnt = mti
    State0%init = mtinit

  end subroutine getStateMT

  ! Set the current state vector and counter
  subroutine setStateMT( State0 )

    implicit none
    type(StateMT_t), intent( in )   :: State0

    if (State0%cnt < 1_wi .or. State0%cnt > n + 1_wi) then
       print *,'err[RandomMT::setStateMT]: state%cnt= ',State0%cnt
       print *,'Must be between 1 and ',n+1_wi
       stop
    endif

    mt = State0%vec
    mti = State0%cnt
    mtinit = State0%init

    if ( .not. State0%init ) then
       print *,'Warning[RandomMT::setStateMT]: mtinit has been set'
       print *,'to false, resetting state%cnt to ',n+1_wi
       mti = n + 1_wi       
    endif

  end subroutine setStateMT

  subroutine nextRandMT(state,ranval)

!   A repeatable sequence of random numbers can be obtained by initializing
!   the state with a call to initRandMT and then getting the state with
!   getRandMT.  This sequence repeats independently of any calls to nextRandMT 
!   or randMT that may occur with other states.

    type(StateMT_t), intent( inout )   :: state
    real, intent( inout )              :: ranval

    if ( .not. state%init ) then
      print *,'err[RandomMT::nexRandMT]: state has not been initialized.'
      print *,'Check whether state was initialized and returned with '
      print *,'initRandMT(iseed,state).'
      stop
    endif

    call setStateMT(state)
    ranval = randMT()
    call getStateMT(state)
    
  end subroutine nextRandMT

  subroutine nextRandMTd(state,ranval)

!   A repeatable sequence of random numbers can be obtained by initializing
!   the state with a call to initRandMT and then getting the state with
!   getRandMT.  This sequence repeats independently of any calls to nextRandMT 
!   or randMT that may occur with other states.

    type(StateMT_t), intent( inout )   :: state
    double precision, intent( inout )  :: ranval

    if ( .not. state%init ) then
      print *,'err[RandomMT::nexRandMT]: state has not been initialized.'
      print *,'Check whether state was initialized and returned with '
      print *,'initRandMT(iseed,state).'
      stop
    endif

    call setStateMT(state)
    ranval =  randMTd()
    call getStateMT(state)
    
  end subroutine nextRandMTd

end module RandomMT
