      subroutine idlmvltlnd(argc,argv)
c     IDL call must be lower case ('idlconstinit_') even though upper is
c     used here

      Implicit None
      Integer*4 argc,argv(*),j
       	
      j = LOC(argc)
      
      call idlmvltlnd_core(%val(argv(1)),
     +     %val(argv(2)),
     +     %val(argv(3)),
     +     %val(argv(4)),
     +     %val(argv(5)),
     +     %val(argv(6)))
      return
      end

c*****************************************************************************
      subroutine idlmvltlnd_core(rlat1,elon1,dismv,dirmv,rlat2,elon2)
      
      use GeomSpheric

      implicit none
      
      real :: rlat1,elon1,dismv,dirmv,rlat2,elon2
      real :: rlon1,rlon2,dirmvc

! convert input direction from 0=N, +clockwise to 0=S, +ccw
      dirmvc=180.-dirmv
      if (dirmvc > 180.)dirmvc=dirmvc-360.
      if (dirmvc < -180.)dirmvc=dirmvc+360.
! convert input lon from east to west
      rlon1=-elon1
      call mvltlnd(rlat1,rlon1,dismv,dirmvc,rlat2,rlon2)
! convert output lon from west to east
      elon2=-rlon2
      
      return
      end
