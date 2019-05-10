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

!--------------------------------------------------------------
!
!  MODULE Invert: contains subroutines:
!                 invrt2_cholesky, invrt1_cholesky, invrt1_LM_cholesky, 
!                 invrt1_nonDiagSy_cholesky, invrt1_LM_nonDiagSy_cholesky
!                 invrt2_GaussJordan, invrt1_GaussJordan, dminv and dmxd
!  Note: The 4 variations of invrt1 with Cholesky decomposition are currently
!        separate codes, although they are nearly identical.  
!        This is an interim condition until they are merged (with optional
!        arguments) and timing tests determine no significant impact.
!        For now, they must be kept consistent manually.
!
!     Copyright AER, Inc., 2002, 2003. All rights Reserved.
!--------------------------------------------------------------
MODULE InvertModule

! <f90Module>***********************************************************
!
! NAME:
!
!   InvertModule
!
! PURPOSE:
!
!   Contains inversion subroutines, retrieving a vector of variables from a 
!   vector of measurements.
!
! INCLUDES:
!
!   None
!
!***********************************************************</f90Module>

  IMPLICIT NONE

  PRIVATE
  PUBLIC  :: invrt1, invrt2
  PUBLIC  :: invrt1_cholesky, invrt2_cholesky
  PUBLIC  :: invrt1_LM_cholesky, invrt1_nonDiagSy_cholesky
  PUBLIC  :: invrt1_LM_nonDiagSy_cholesky
  PUBLIC  :: invrt1_GaussJordan, invrt2_GaussJordan
  PUBLIC  :: invrt1_GJ_factors
  PUBLIC  :: dminv, invrtRegr

  interface invrt1
    module procedure invrt1_cholesky
  end interface

  interface invrt1_LM
    module procedure invrt1_LM_cholesky
  end interface

  interface invrt1_LM_nonDiagSy
    module procedure invrt1_LM_nonDiagSy_cholesky
  end interface

  interface invrt1_nonDiagSy
    module procedure invrt1_nonDiagSy_cholesky
  end interface

  interface invrt2
    module procedure invrt2_cholesky
  end interface


  !---------------------------------------------------------------
  !     Declarations of private data (private to the module)
  !---------------------------------------------------------------
  INTEGER, parameter :: DBL=SELECTED_REAL_KIND(15)

  INTEGER, PUBLIC, parameter :: DBL_INVT=DBL   ! For consistency of arguments 
                                               ! in external calls 
        
CONTAINS

  SUBROUTINE invrt2_cholesky(xkt,e,bcov,ydif,xn1,xn,np,nchan)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrt2_cholesky
!
! PURPOSE:
!
!   Performs profile inversion by optimal estimation in configuration 2, using 
!   Cholesky decomposition.
!
! SYNTAX:
!
!   CALL invrt2_cholesky(xkt, e, bcov, ydif, xn1, xn, np, nchan)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xkt    REAL     Jacobian matrix transpose 
!   e      REAL     measurement error covariance matrix
!   bcov   REAL     Background error covariance matrix
!   ydif   REAL     Radiometric measurements minus computed
!   xn1    REAL     initial profile
!   np     INTEGER  number of parameters
!   nchan  INTEGER  Number of channels
!   
!   INPUTS/OUTPUTS:
!   
!   xn     REAL     inverted profile
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  ! Function name: invrt2_cholesky
  ! Purpose: Performs profile inversion.
  ! Description: The inversion is performed by non-linear optimal
  !              estimation using a forecast profile and its error
  !              covariance as backgroung constraints. Please refer
  !              to "CrIS ATBD" by AER in detail.
  !              In order to accomodate the cray convention 
  !              (row dimension = leading array dimension), 
  !              the k matrix passed as an argument is stored
  !              as the transpose of the k matrix in xkt. The other 
  !              matrices (e,bcov) are symmetrical. This version of invrt 
  !              is to be used when the number of channels is small 
  !              compared to the number of parameters.
  ! Inputs:
  ! Var_Name      Type         Description
  ! --------      ----         -----------
  ! xkt           real         matrix to be inverted.
  ! e             real         error covariance matrix.
  ! bcov          real         profile covariance matrix.
  ! ydif          real         difference of radiances.
  ! xn1           real         initial profile.
  ! np            real         number of parameters.
  ! nchan         real         number of channels.     
  !
  ! Outputs:
  ! Var_Name      Type       Description
  ! --------      ----       -----------
  ! xn            real       inverted profile
  ! SxMwOut       real       covariance matrix in EOF domain
  !                          (used for cascade mode)  
  !
  ! Copyright: Atmospheric and Environmental Research, Inc., 1997        
  ! Developed by Atmospheric and Environmental Research, Inc.            
  ! Modified, March, 2003 by W. Gallery, AER, Inc.
  !    Adapted to F90 syntax
  !    Replaced matrix subroutines with intrinsics
  ! Modified April, 2006 by V. Payne, AER, Inc
  !    Implemented Cholesky decomposition instead of dminv
  !--------------------------------------------------------------
 
    !----Input Variables
    INTEGER,                    INTENT(IN)    :: np, nchan
    REAL,  DIMENSION(np,nchan), INTENT(IN)    :: xkt
    REAL,  DIMENSION(nchan),    INTENT(IN)    :: e,ydif
    REAL,  DIMENSION(np,np),    INTENT(IN)    :: bcov
    REAL,  DIMENSION(np),       INTENT(IN)    :: xn1
    REAL,  DIMENSION(np),       INTENT(INOUT) :: xn

    !----Variables of user-defined precision
    REAL(KIND=DBL), DIMENSION(np,nchan)     :: dxkt,dwk
    REAL(KIND=DBL), DIMENSION(nchan)        :: de,dydif,dyx
    REAL(KIND=DBL), DIMENSION(np,np)        :: dc
    REAL(KIND=DBL), DIMENSION(np)           :: dxn1,dxn
    REAL(KIND=DBL), DIMENSION(nchan,np)     :: dxk
    REAL(KIND=DBL), DIMENSION(nchan,nchan)  :: dh
    REAL(KIND=DBL), DIMENSION(nchan)        :: z

    INTEGER                                 :: i
    INTEGER                                 :: info, info2

    !----Conversion from REAL to user-defined REAL(KIND=DBL)
    dxkt=REAL(xkt,KIND=DBL)
    dc=REAL(bcov,KIND=DBL)
    de=REAL(e,KIND=DBL)
    dxn1=REAL(xn1,KIND=DBL)
    dydif=REAL(ydif,KIND=DBL)

    !----Transpose k matrix
    dxk=TRANSPOSE(dxkt) 
    !----Contribution dy=h-1ktw as g-1kt(kg-1kt+w-1)-1.
    dwk=MATMUL(dc, dxkt)
    dh=MATMUL(dxk, dwk)
    !----Replace diagonal of dh with sum of diagonal of dh & de
    DO i=1,nchan
       dh(i,i)=dh(i,i)+de(i)
    END DO

    dyx = dydif + MATMUL(dxk, dxn1)
    
    !duplicate dyx array for input/output to Cholesky routine
    z = dyx
     
    ! LAPACK routine for Cholesky decomposition
    CALL dpotrf('L',nchan,dh,nchan,info)
    IF (info .NE. 0) THEN
       print*, 'err[Invert::invrt2_cholesky]:'
       print*, 'Problem with Cholesky decomposition.'
       CALL errorHalt(1)
    ENDIF
    CALL dpotrs('L',nchan,1,dh,nchan,z,nchan,info2)
    IF (info2 .NE. 0) THEN
       print*, 'err[Invert::invrt2_cholesky]:'
       print*, 'Problem with solution of inverse.'
       CALL errorHalt(1)
    ENDIF
    
    !---update to state vector
    dxn = MATMUL(dwk,z)

    
    !----Conversion from REAL(KIND=DBL) to REAL
    xn=REAL(dxn)

    RETURN
  END SUBROUTINE invrt2_cholesky

  subroutine invrt1_cholesky(xkt,e,cinv,ydif,xn1,xn,np,nchan,hm1,dofs)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrt1_cholesky
!
! PURPOSE:
!
!   Performs profile inversion by optimal estimation in configuration 1, using 
!   Cholesky decomposition.
!
! SYNTAX:
!
!   CALL invrt1_cholesky(xkt, e, cinv, ydif, xn1, xn, np, nchan)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xkt    REAL     Jacobian matrix transpose 
!   e      REAL     measurement error covariance matrix
!   cinv   REAL     inverse of background error covariance matrix
!   ydif   REAL     Radiometric measurements minus computed
!   xn1    REAL     initial profile
!   np     INTEGER  number of parameters
!   nchan  INTEGER  Number of channels
!   
!   INPUTS/OUTPUTS:
!   
!   xn     REAL     inverted profile
!   hm1    REAL     retrieved error covariance in retrieval space
!   dofs   REAL     degrees of freedom for signal
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !***********************************************************************
    !* Function name: invrt1_cholesky
    !* Purpose: Performs profile inversion
    !* Description: The inversion is performed by non-linear optimal
    !*              estimation using a background profile and its error
    !*              covariance as constraints (please refer to the ATBD 
    !*              for details). In order to accomodate the cray 
    !*              convention (row dimension = leading array dimension), 
    !*              the k matrix passed as an argument is stored
    !*              as the transpose of the k matrix in xkt. The 
    !*              other matrices (e,cinv) are symmetrical. This version 
    !*              of invrt is to be used when the number of channels 
    !*              is larger than the number of retrieved parameters
    !*              (i.e., in the MW+IR retrieval). All vector and matrix
    !*              operations are in double precision.
    !* Inputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xkt           real       matrix to be inverted.
    !* e             real       measurement error covariance matrix.
    !* cinv          real       inverse of background profile covariance matrix.
    !* ydif          real       difference of radiances.
    !* xn1           real       initial (guess) profile.
    !* np            real       number of retrieved parameters.
    !* nchan         real       number of channels.  
    !*
    !* Outputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xn            real       inverted profile
    !* h             real       inverse of retrieved error covariance in retrieval space
    !*
    !* Externals: 
    !* Name    Description
    !* ----    -----------
    !* mxd     multiplication of a matrix with each element of a vector
    !* dminv   matrix inversion
    !*
    !*<copyright>=========================================================
    !*
    !* Developed by Atmospheric and Environmental Research, Inc.                
    !*
    !* Copyright: Atmospheric and Environmental Inc., 1997-2004
    !*                 All Rights Reserved
    !========================================================</copyright>
    !***********************************************************************

    !----Input Variables
    integer,                    intent(in)    :: np, nchan
    real,  dimension(np,nchan), intent(in)    :: xkt
    real,  dimension(nchan),    intent(in)    :: e,ydif
    real,  dimension(np,np),    intent(in)    :: cinv
    real,  dimension(np),       intent(in)    :: xn1
    real,  dimension(np),       intent(inout) :: xn
    real,  dimension(np,np),    intent(inout) :: hm1
    real,                       intent(inout) :: dofs
    
    !----Variables of user-defined precision
    real(kind=DBL), dimension(np,nchan)     :: dxkt,dwk
    real(kind=DBL), dimension(nchan)        :: de,dydif,dyx
    real(kind=DBL), dimension(np,np)        :: dcinv,dhm1,dak
    real(kind=DBL), dimension(np)           :: dxn1,dxn,dkyx,ddiag
    real(kind=DBL), dimension(nchan,np)     :: dxk
    real(kind=DBL), dimension(np,np)        :: dh,dktsmk
    real(kind=DBL)                          :: ddet
    
    INTEGER                                 :: info,info2,ip

    !----Conversion from real to user-defined real(kind=DBL)
    dxkt=real(xkt,kind=DBL)
    dcinv=real(cinv,kind=DBL)
    de=real(1./e,kind=DBL)
    dxn1=real(xn1,kind=DBL)
    dydif=real(ydif,kind=DBL)

    print *,'ydif ',ydif
    do ip=1,10
       print *,ip,xn1(ip)
    end do
    stop 

    !----Transpose k matrix
    dxk=transpose(dxkt) 
    !--   calculate contribution functions dy=h-1ktw as (ktwk+g)-1ktw
    call dmxd(dxkt,de,dwk,np,nchan)
    dktsmk = matmul(dwk,dxk)
    dh = dktsmk+dcinv

! New code for error covariance and AK
    dhm1=dh
    CALL dminv(dhm1,np,ddet)

    dak = matmul(dhm1,dktsmk)
    do ip=1,np
       ddiag(ip) = dak(ip,ip)
    enddo

    dyx = dydif + MATMUL(dxk,dxn1)
    dkyx = MATMUL(dwk,dyx)
    dxn = dkyx

    CALL dpotrf('L',np,dh,np,info)
    IF (info .NE. 0) THEN
       print*, 'err[Invert::invrt1_cholesky]:'
       print*, 'Problem with Cholesky decomposition.'
       CALL errorHalt(1)
    ENDIF
    CALL dpotrs('L',np,1,dh,np,dxn,np,info2)
    IF (info2 .NE. 0) THEN
       print*, 'err[Invert::invrt1_cholesky]:'
       print*, 'Problem with solution of inverse.'
       CALL errorHalt(1)
    ENDIF

    !----Conversion from real(kind=DBL) to real
    xn=real(dxn)
    hm1 = real(dhm1)
    dofs = real(sum(ddiag))

    return
  end subroutine invrt1_cholesky

  subroutine invrt1_LM_cholesky(xkt,e,cinv,ydif,xn1,xn,np,nchan,gm)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrt1_LM_cholesky
!
! PURPOSE:
!
!   Performs profile inversion by optimal estimation in configuration 1, using 
!   Cholesky decomposition, with Levenberg-Marquardt
!
! SYNTAX:
!
!   CALL invrt1_LM_cholesky(xkt, e, cinv, ydif, xn1, xn, np, nchan, gm)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xkt    REAL     Jacobian matrix transpose 
!   e      REAL     measurement error covariance matrix
!   cinv   REAL     inverse of background error covariance matrix
!   ydif   REAL     Radiometric measurements minus computed
!   xn1    REAL     initial profile
!   np     INTEGER  number of parameters
!   nchan  INTEGER  Number of channels
!   gm     REAL     Levenberg-Marquardt parameter
!   
!   INPUTS/OUTPUTS:
!   
!   xn     REAL     inverted profile
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !***********************************************************************
    !* Function name: invrt1_LM_cholesky
    !* Purpose: Performs profile inversion with Levenberg-Marquardt
    !* Description: The inversion is performed by non-linear optimal
    !*              estimation using a background profile and its error
    !*              covariance as constraints (please refer to the ATBD 
    !*              for details). In order to accomodate the cray 
    !*              convention (row dimension = leading array dimension), 
    !*              the k matrix passed as an argument is stored
    !*              as the transpose of the k matrix in xkt. The 
    !*              other matrices (e,cinv) are symmetrical. This version 
    !*              of invrt is to be used when the number of channels 
    !*              is larger than the number of retrieved parameters
    !*              (i.e., in the MW+IR retrieval). All vector and matrix
    !*              operations are in double precision.
    !* Inputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xkt           real       matrix to be inverted.
    !* e             real       measurement error covariance matrix.
    !* cinv          real       inverse of background profile covariance matrix.
    !* ydif          real       difference of radiances.
    !* xn1           real       initial (guess) profile.
    !* np            real       number of retrieved parameters.
    !* nchan         real       number of channels.  
    !* gm            real       Levenberg-Marquardt parameter
    !*
    !* Outputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xn            real       inverted profile
    !*
    !* Externals: 
    !* Name    Description
    !* ----    -----------
    !* mxd     multiplication of a matrix with each element of a vector
    !* dminv   matrix inversion
    !*
    !*<copyright>=========================================================
    !*
    !* Developed by Atmospheric and Environmental Research, Inc.                
    !*
    !* Copyright: Atmospheric and Environmental Inc., 1997-2004
    !*                 All Rights Reserved
    !========================================================</copyright>
    !***********************************************************************

    !----Input Variables
    integer,                    intent(in)    :: np, nchan
    real,  dimension(np,nchan), intent(in)    :: xkt
    real,  dimension(nchan),    intent(in)    :: e,ydif
    real,  dimension(np,np),    intent(in)    :: cinv
    real,  dimension(np),       intent(in)    :: xn1
    real,  dimension(np),       intent(inout) :: xn
    real,                       intent(in)    :: gm
    
    !----Variables of user-defined precision
    real(kind=DBL), dimension(np,nchan)     :: dxkt,dwk
    real(kind=DBL), dimension(nchan)        :: de,dydif,dyx
    real(kind=DBL), dimension(np,np)        :: dcinv
    real(kind=DBL), dimension(np)           :: dxn1,dxn,dkyx
    real(kind=DBL), dimension(nchan,np)     :: dxk
    real(kind=DBL), dimension(np,np)        :: dh
    
    INTEGER                                 :: info,info2

    !----Conversion from real to user-defined real(kind=DBL)
    dxkt=real(xkt,kind=DBL)
    dcinv=real(cinv,kind=DBL)
    de=real(1./e,kind=DBL)
    dxn1=real(xn1,kind=DBL)
    dydif=real(ydif,kind=DBL)

    !----Transpose k matrix
    dxk=transpose(dxkt) 
    !--   calculate contribution functions dy=h-1ktw as (ktwk+g)-1ktw
    call dmxd(dxkt,de,dwk,np,nchan)
    dh = matmul(dwk,dxk)+dcinv*(1+gm)

    dyx = dydif + 1./(1+gm)*MATMUL(dxk,dxn1)

    dkyx = MATMUL(dwk,dyx)

    dxn = dkyx

    CALL dpotrf('L',np,dh,np,info)
    IF (info .NE. 0) THEN
       print*, 'err[Invert::invrt1_LM_cholesky]:'
       print*, 'Problem with Cholesky decomposition.'
       CALL errorHalt(1)
    ENDIF
    CALL dpotrs('L',np,1,dh,np,dxn,np,info2)
    IF (info2 .NE. 0) THEN
       print*, 'err[Invert::invrt1_LM_cholesky]:'
       print*, 'Problem with solution of inverse.'
       CALL errorHalt(1)
    ENDIF

    !----Conversion from real(kind=DBL) to real
    xn=real(dxn+gm/(1.+gm)*dxn1)

    return
  end subroutine invrt1_LM_cholesky

  subroutine invrt1_nonDiagSy_cholesky(xkt,einv,cinv,ydif,xn1,xn,np,nchan)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrt1_nonDiagSy_cholesky
!
! PURPOSE:
!
!   Performs profile inversion by optimal estimation in configuration 1, using 
!   Cholesky decomposition, where measurement error covariance matrix is not 
!   assumed to be diagonal.
!
! SYNTAX:
!
!   CALL invrt1_nonDiagSy_cholesky(xkt, einv, cinv, ydif, xn1, xn, 
!      np, nchan)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xkt    REAL     Jacobian matrix transpose 
!   einv   REAL     inverse of measurement error covariance matrix.
!   cinv   REAL     inverse of background error covariance matrix
!   ydif   REAL     Radiometric measurements minus computed
!   xn1    REAL     initial profile
!   np     INTEGER  number of parameters
!   nchan  INTEGER  Number of channels
!   
!   INPUTS/OUTPUTS:
!   
!   xn     REAL     inverted profile
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !***********************************************************************
    !* Function name: invrt1_nonDiagSy_cholesky
    !* Purpose: Performs profile inversion, where
    !*          measurement error covariance matrix is not assumed
    !*          to be diagonal.
    !* Description: The inversion is performed by non-linear optimal
    !*              estimation using a background profile and its error
    !*              covariance as constraints (please refer to the ATBD 
    !*              for details). In order to accomodate the cray 
    !*              convention (row dimension = leading array dimension), 
    !*              the k matrix passed as an argument is stored
    !*              as the transpose of the k matrix in xkt. The 
    !*              other matrices (einv,cinv) are symmetrical. This version 
    !*              of invrt is to be used when the number of channels 
    !*              is larger than the number of retrieved parameters
    !*              (i.e., in the MW+IR retrieval). All vector and matrix
    !*              operations are in double precision.
    !* Inputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xkt           real       matrix to be inverted.
    !* einv          real       inverse of measurement error covariance matrix.
    !* cinv          real       inverse of background profile covariance matrix.
    !* ydif          real       difference of radiances.
    !* xn1           real       initial (guess) profile.
    !* np            real       number of retrieved parameters.
    !* nchan         real       number of channels.  
    !*
    !* Outputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xn            real       inverted profile
    !*
    !* Externals: 
    !* Name    Description
    !* ----    -----------
    !* mxd     multiplication of a matrix with each element of a vector
    !* dminv   matrix inversion
    !*
    !*<copyright>=========================================================
    !*
    !* Developed by Atmospheric and Environmental Research, Inc.                
    !*
    !* Copyright: Atmospheric and Environmental Inc., 1997-2004
    !*                 All Rights Reserved
    !========================================================</copyright>
    !***********************************************************************

    !----Input Variables
    integer,                       intent(in)    :: np, nchan
    real,  dimension(np,nchan),    intent(in)    :: xkt
    real,  dimension(nchan),       intent(in)    :: ydif
    real,  dimension(nchan,nchan), intent(in)    :: einv
    real,  dimension(np,np),       intent(in)    :: cinv
    real,  dimension(np),          intent(in)    :: xn1
    real,  dimension(np),          intent(inout) :: xn
    
    !----Variables of user-defined precision
    real(kind=DBL), dimension(np,nchan)     :: dxkt,dwk
    real(kind=DBL), dimension(nchan)        :: dydif,dyx
    real(kind=DBL), dimension(nchan,nchan)  :: deinv
    real(kind=DBL), dimension(np,np)        :: dcinv
    real(kind=DBL), dimension(np)           :: dxn1,dxn,dkyx
    real(kind=DBL), dimension(nchan,np)     :: dxk
    real(kind=DBL), dimension(np,np)        :: dh
    
    INTEGER                                 :: info,info2

    !----Conversion from real to user-defined real(kind=DBL)
    dxkt=real(xkt,kind=DBL)
    dcinv=real(cinv,kind=DBL)
    deinv=real(einv,kind=DBL)
    dxn1=real(xn1,kind=DBL)
    dydif=real(ydif,kind=DBL)

    !----Transpose k matrix
    dxk=transpose(dxkt) 
    !--   calculate contribution functions dy=h-1ktw as (ktwk+g)-1ktw
    dwk=matmul(dxkt,deinv)
    dh = matmul(dwk,dxk)+dcinv

    dyx = dydif + MATMUL(dxk,dxn1)

    dkyx = MATMUL(dwk,dyx)

    dxn = dkyx

    CALL dpotrf('L',np,dh,np,info)
    IF (info .NE. 0) THEN
       print*, 'err[Invert::invrt1_nonDiagSy_cholesky]:'
       print*, 'Problem with Cholesky decomposition.'
       CALL errorHalt(1)
    ENDIF
    CALL dpotrs('L',np,1,dh,np,dxn,np,info2)
    IF (info2 .NE. 0) THEN
       print*, 'err[Invert::invrt1_nonDiagSy_cholesky]:'
       print*, 'Problem with solution of inverse.'
       CALL errorHalt(1)
    ENDIF

    !----Conversion from real(kind=DBL) to real
    xn=real(dxn)

    return
  end subroutine invrt1_nonDiagSy_cholesky

  subroutine invrt1_LM_nonDiagSy_cholesky(xkt,einv,cinv,ydif,xn1,xn,np,nchan,gm)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrt1_LM_nonDiagSy_cholesky
!
! PURPOSE:
!
!   Performs profile inversion by optimal estimation in configuration 1, using 
!   Cholesky decomposition, with Levenberg-Marquardt, and measurement error 
!   covariance matrix is not assumed to be diagonal.
!
! SYNTAX:
!
!   CALL invrt1_LM_nonDiagSy_cholesky(xkt, einv, cinv, ydif, xn1, xn, 
!      np, nchan, gm)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xkt    REAL     Jacobian matrix transpose 
!   einv   REAL     inverse of measurement error covariance matrix.
!   cinv   REAL     inverse of background error covariance matrix
!   ydif   REAL     Radiometric measurements minus computed
!   xn1    REAL     initial profile
!   np     INTEGER  number of parameters
!   nchan  INTEGER  Number of channels
!   gm     REAL     Levenberg-Marquardt parameter
!   
!   INPUTS/OUTPUTS:
!   
!   xn     REAL     inverted profile
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !***********************************************************************
    !* Function name: invrt1_LM_nonDiagSy_cholesky
    !* Purpose: Performs profile inversion with Levenberg-Marquardt, 
    !*          and measurement error covariance matrix is not assumed
    !*          to be diagonal.
    !* Description: The inversion is performed by non-linear optimal
    !*              estimation using a background profile and its error
    !*              covariance as constraints (please refer to the ATBD 
    !*              for details). In order to accomodate the cray 
    !*              convention (row dimension = leading array dimension), 
    !*              the k matrix passed as an argument is stored
    !*              as the transpose of the k matrix in xkt. The 
    !*              other matrices (einv,cinv) are symmetrical. This version 
    !*              of invrt is to be used when the number of channels 
    !*              is larger than the number of retrieved parameters
    !*              (i.e., in the MW+IR retrieval). All vector and matrix
    !*              operations are in double precision.
    !* Inputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xkt           real       matrix to be inverted.
    !* einv          real       inverse of measurement error covariance matrix.
    !* cinv          real       inverse of background profile covariance matrix.
    !* ydif          real       difference of radiances.
    !* xn1           real       initial (guess) profile.
    !* np            real       number of retrieved parameters.
    !* nchan         real       number of channels.  
    !* gm            real       Levenberg-Marquardt parameter
    !*
    !* Outputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xn            real       inverted profile
    !*
    !* Externals: 
    !* Name    Description
    !* ----    -----------
    !* mxd     multiplication of a matrix with each element of a vector
    !* dminv   matrix inversion
    !*
    !*<copyright>=========================================================
    !*
    !* Developed by Atmospheric and Environmental Research, Inc.                
    !*
    !* Copyright: Atmospheric and Environmental Inc., 1997-2004
    !*                 All Rights Reserved
    !========================================================</copyright>
    !***********************************************************************

    !----Input Variables
    integer,                       intent(in)    :: np, nchan
    real,  dimension(np,nchan),    intent(in)    :: xkt
    real,  dimension(nchan),       intent(in)    :: ydif
    real,  dimension(nchan,nchan), intent(in)    :: einv
    real,  dimension(np,np),       intent(in)    :: cinv
    real,  dimension(np),          intent(in)    :: xn1
    real,  dimension(np),          intent(inout) :: xn
    real,                          intent(in)    :: gm
    
    !----Variables of user-defined precision
    real(kind=DBL), dimension(np,nchan)     :: dxkt,dwk
    real(kind=DBL), dimension(nchan)        :: dydif,dyx
    real(kind=DBL), dimension(nchan,nchan)  :: deinv
    real(kind=DBL), dimension(np,np)        :: dcinv
    real(kind=DBL), dimension(np)           :: dxn1,dxn,dkyx
    real(kind=DBL), dimension(nchan,np)     :: dxk
    real(kind=DBL), dimension(np,np)        :: dh
    
    INTEGER                                 :: info,info2

    !----Conversion from real to user-defined real(kind=DBL)
    dxkt=real(xkt,kind=DBL)
    dcinv=real(cinv,kind=DBL)
    deinv=real(einv,kind=DBL)
    dxn1=real(xn1,kind=DBL)
    dydif=real(ydif,kind=DBL)

    !----Transpose k matrix
    dxk=transpose(dxkt) 
    !--   calculate contribution functions dy=h-1ktw as (ktwk+g)-1ktw
    dwk=matmul(dxkt,deinv)
    dh = matmul(dwk,dxk)+dcinv*(1+gm)

    dyx = dydif + 1./(1+gm)*MATMUL(dxk,dxn1)

    dkyx = MATMUL(dwk,dyx)

    dxn = dkyx

    CALL dpotrf('L',np,dh,np,info)
    IF (info .NE. 0) THEN
       print*, 'err[Invert::invrt1_LM_nonDiagSy_cholesky]:'
       print*, 'Problem with Cholesky decomposition.'
       CALL errorHalt(1)
    ENDIF
    CALL dpotrs('L',np,1,dh,np,dxn,np,info2)
    IF (info2 .NE. 0) THEN
       print*, 'err[Invert::invrt1_LM_nonDiagSy_cholesky]:'
       print*, 'Problem with solution of inverse.'
       CALL errorHalt(1)
    ENDIF

    !----Conversion from real(kind=DBL) to real
    xn=real(dxn+gm/(1.+gm)*dxn1)

    return
  end subroutine invrt1_LM_nonDiagSy_cholesky

  SUBROUTINE invrt2_GaussJordan(xkt,e,bcov,ydif,xn1,xn,np,nchan)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrt2_GaussJordan
!
! PURPOSE:
!
!   Performs profile inversion by optimal estimation in configuration 2, using 
!   Gaus-Jordan method.
!
! SYNTAX:
!
!   CALL invrt2_GaussJordan(xkt, e, bcov, ydif, xn1, xn, np, nchan)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xkt    REAL     Jacobian matrix transpose 
!   e      REAL     measurement error covariance matrix
!   bcov   REAL     Background error covariance matrix
!   ydif   REAL     Radiometric measurements minus computed
!   xn1    REAL     initial profile
!   np     INTEGER  number of parameters
!   nchan  INTEGER  Number of channels
!   
!   INPUTS/OUTPUTS:
!   
!   xn     REAL     inverted profile
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  ! Function name: invrt2_GaussJordan
  ! Purpose: Performs profile inversion.
  ! Usage: call invrt2_GaussJordan(xkt,e,bcov,ydif,xn1,xn,np,nchan)
  ! Description: The inversion is performed by non-linear optimal
  !              estimation using a forecast profile and its error
  !              covariance as backgroung constraints. Please refer
  !              to "CrIS ATBD" by AER in detail.
  !              In order to accomodate the cray convention 
  !              (row dimension = leading array dimension), 
  !              the k matrix passed as an argument is stored
  !              as the transpose of the k matrix in xkt. The other 
  !              matrices (e,bcov) are symmetrical. This version of invrt 
  !              is to be used when the number of channels is small 
  !              compared to the number of parameters.
  ! Inputs:
  ! Var_Name      Type         Description
  ! --------      ----         -----------
  ! xkt           real         matrix to be inverted.
  ! e             real         error covariance matrix.
  ! bcov          real         profile covariance matrix.
  ! ydif          real         difference of radiances.
  ! xn1           real         initial profile.
  ! np            real         number of parameters.
  ! nchan         real         number of channels.     
  !
  ! Outputs:
  ! Var_Name      Type       Description
  ! --------      ----       -----------
  ! xn            real       inverted profile
  ! SxMwOut       real       covariance matrix in EOF domain
  !                          (used for cascade mode)  
  !
  ! Copyright: Atmospheric and Environmental Research, Inc., 1997        
  ! Developed by Atmospheric and Environmental Research, Inc.            
  ! Modified, March, 2003 by W. Gallery, AER, Inc.
  !    Adapted to F90 syntax
  !    Replaced matrix subroutines with intrinsics
  !--------------------------------------------------------------
 
    !----Input Variables
    INTEGER,                    INTENT(IN)    :: np, nchan
    REAL,  DIMENSION(np,nchan), INTENT(IN)    :: xkt
    REAL,  DIMENSION(nchan),    INTENT(IN)    :: e,ydif
    REAL,  DIMENSION(np,np),    INTENT(IN)    :: bcov
    REAL,  DIMENSION(np),       INTENT(IN)    :: xn1
    REAL,  DIMENSION(np),       INTENT(INOUT) :: xn

    !----Variables of user-defined precision
    REAL(KIND=DBL), DIMENSION(np,nchan)     :: dxkt,dwk,ddy
    REAL(KIND=DBL), DIMENSION(nchan)        :: de,dydif
    REAL(KIND=DBL), DIMENSION(np,np)        :: dc,da
    REAL(KIND=DBL), DIMENSION(np)           :: dxn1,dxn
    REAL(KIND=DBL), DIMENSION(nchan,np)     :: dxk
    REAL(KIND=DBL), DIMENSION(nchan,nchan)  :: dh
    REAL(KIND=DBL)                          :: ddet

    INTEGER                                 :: i

    !----Conversion from REAL to user-defined REAL(KIND=DBL)
    dxkt=REAL(xkt,KIND=DBL)
    dc=REAL(bcov,KIND=DBL)
    de=REAL(e,KIND=DBL)
    dxn1=REAL(xn1,KIND=DBL)
    dydif=REAL(ydif,KIND=DBL)

    !----Transpose k matrix
    dxk=TRANSPOSE(dxkt) 
    !----Contribution dy=h-1ktw as g-1kt(kg-1kt+w-1)-1.
    dwk=MATMUL(dc, dxkt)
    dh=MATMUL(dxk, dwk)
    !----Replace diagonal of dh with sum of diagonal of dh & de
    DO i=1,nchan
       dh(i,i)=dh(i,i)+de(i)
    END DO
    CALL dminv(dh,nchan,ddet)
    ddy=MATMUL(dwk,dh)
    !----Calculate averaging kernel a=dy*k
    da=MATMUL(ddy, dxk)
    !---New estimate of profile as xn-x0 = a (xn-1-x0)+ dy (ym-yn-1)
    dxn=MATMUL(da, dxn1)+MATMUL(ddy, dydif)

    !----Conversion from REAL(KIND=DBL) to REAL
    xn=REAL(dxn)

    RETURN
  END SUBROUTINE invrt2_GaussJordan

  subroutine invrt1_GaussJordan(xkt,e,cinv,ydif,xn1,xn,np,nchan)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrt1_GaussJordan
!
! PURPOSE:
!
!   Performs profile inversion by optimal estimation in configuration 1, using 
!   Gaus-Jordan method.
!
! SYNTAX:
!
!   CALL invrt1_GaussJordan(xkt, e, cinv, ydif, xn1, xn, np, nchan)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xkt    REAL     Jacobian matrix transpose 
!   e      REAL     measurement error covariance matrix
!   cinv   REAL     inverse of background error covariance matrix
!   ydif   REAL     Radiometric measurements minus computed
!   xn1    REAL     initial profile
!   np     INTEGER  number of parameters
!   nchan  INTEGER  Number of channels
!   
!   INPUTS/OUTPUTS:
!   
!   xn     REAL     inverted profile
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !***********************************************************************
    !* Function name: invrt1_GaussJordan
    !* Purpose: Performs profile inversion
    !* Usage: call invrt1_GaussJordan(xkt,e,cinv,ydif,xn1,xn,np,nch)
    !* Description: The inversion is performed by non-linear optimal
    !*              estimation using a background profile and its error
    !*              covariance as constraints (please refer to the ATBD 
    !*              for details). In order to accomodate the cray 
    !*              convention (row dimension = leading array dimension), 
    !*              the k matrix passed as an argument is stored
    !*              as the transpose of the k matrix in xkt. The 
    !*              other matrices (e,cinv) are symmetrical. This version 
    !*              of invrt is to be used when the number of channels 
    !*              is larger than the number of retrieved parameters
    !*              (i.e., in the MW+IR retrieval). All vector and matrix
    !*              operations are in double precision.
    !* Inputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xkt           real       matrix to be inverted.
    !* e             real       measurement error covariance matrix.
    !* cinv          real       inverted background profile covariance matrix.
    !* ydif          real       difference of radiances.
    !* xn1           real       initial (guess) profile.
    !* np            real       number of retrieved parameters.
    !* nchan         real       number of channels.  
    !*
    !* Outputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xn            real       inverted profile
    !* hR            real       retrieval covariance in retrieval space
    !*
    !* Externals: 
    !* Name    Description
    !* ----    -----------
    !* mxd     multiplication of a matrix with each element of a vector
    !* dminv   matrix inversion
    !*
    !*<copyright>=========================================================
    !*
    !* Developed by Atmospheric and Environmental Research, Inc.                
    !*
    !* Copyright: Atmospheric and Environmental Inc., 1997-2004
    !*                 All Rights Reserved
    !========================================================</copyright>
    !***********************************************************************

    !----Input Variables
    integer,                    intent(in)    :: np, nchan
    real,  dimension(np,nchan), intent(in)    :: xkt
    real,  dimension(nchan),    intent(in)    :: e,ydif
    real,  dimension(np,np),    intent(in)    :: cinv
    real,  dimension(np),       intent(in)    :: xn1
    real,  dimension(np),       intent(inout) :: xn
    real,  dimension(np,np)                   :: hR
    
    !----Variables of user-defined precision
    real(kind=DBL), dimension(np,nchan)     :: dxkt,dwk,ddy
    real(kind=DBL), dimension(nchan)        :: de,dydif
    real(kind=DBL), dimension(np,np)        :: dcinv,da
    real(kind=DBL), dimension(np)           :: dxn1,dxn
    real(kind=DBL), dimension(nchan,np)     :: dxk
    real(kind=DBL), dimension(np,np)        :: dh
    real(kind=DBL)                          :: ddet
    
    !----Conversion from real to user-defined real(kind=DBL)
    dxkt=real(xkt,kind=DBL)
    dcinv=real(cinv,kind=DBL)
    de=real(1./e,kind=DBL)
    dxn1=real(xn1,kind=DBL)
    dydif=real(ydif,kind=DBL)

    !----Transpose k matrix
    dxk=transpose(dxkt) 
    !--   calculate contribution functions dy=h-1ktw as (ktwk+g)-1ktw
    call dmxd(dxkt,de,dwk,np,nchan)
    dh = matmul(dwk,dxk)+dcinv
    call dminv(dh,np,ddet)
    ddy=matmul(dh,dwk)
    !----Calculate averaging kernel a=dy*k
    da=matmul(ddy, dxk)
    !---New estimate of profile as xn-x0 = a (xn-1-x0)+ dy (ym-yn-1)
    dxn=matmul(da, dxn1)+matmul(ddy, dydif)

    !----Conversion from real(kind=DBL) to real
    xn=real(dxn)
    !-- Fill in output cov matrix
    hR=real(dh)
    
    return
  end subroutine invrt1_GaussJordan

  subroutine invrt1_GJ_factors(xkt,e,cinv,yges,xn1,np,nchan,xOff,xGain)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrt1_GJ_factors
!
! PURPOSE:
!
!   Compute factors that can be used in inversion by optimal estimation in 
!   configuration 1, using Gaus-Jordan method.
!
! SYNTAX:
!
!   CALL invrt1_GJ_factors(xkt, e, cinv, yges, xn1, np, nchan, xOff, xGain)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   xkt    REAL     Jacobian matrix transpose 
!   e      REAL     measurement error covariance matrix
!   cinv   REAL     inverse of background error covariance matrix
!   yges   REAL     Radiometric data computed from state vector
!   xn1    REAL     initial profile
!   np     INTEGER  number of parameters
!   nchan  INTEGER  Number of channels
!   
!   INPUTS/OUTPUTS:
!   
!   xOff   REAL     linear offset coefficients
!   xGain  REAL     linear gain coefficients
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !***********************************************************************
    !* Function name: invrt1_GJ_factors
    !* Purpose: Performs profile inversion
    !* Usage: call invrt1_GJ_factors(xkt,e,cinv,ydif,xn1,np,nchan,xOff,xGain)
    !* Description: Compute factors that can be used in inversion.
    !*              This is equivalent to invrt1_GaussJordan, except that
    !*              1) it does not complete the inversion, 2) it processes
    !*              the computed radiance from the prior guess (Y) 
    !*              rather than the difference Ym-Y, where Ym is the 
    !*              measured radiance, and 3) the inversion can be
    !*              completed by subsequently combining Ym with output 
    !*              xOff and xGain as x=xOff+xGain*Y.
    !*              The inversion is performed by non-linear optimal
    !*              estimation using a background profile and its error
    !*              covariance as constraints (please refer to the ATBD 
    !*              for details). In order to accomodate the cray 
    !*              convention (row dimension = leading array dimension), 
    !*              the k matrix passed as an argument is stored
    !*              as the transpose of the k matrix in xkt. The 
    !*              other matrices (e,cinv) are symmetrical. This version 
    !*              of invrt is to be used when the number of channels 
    !*              is larger than the number of retrieved parameters
    !*              (i.e., in the MW+IR retrieval). All vector and matrix
    !*              operations are in double precision.
    !* Inputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xkt           real       matrix to be inverted.
    !* e             real       measurement error covariance matrix.
    !* cinv          real       inverted background profile covariance matrix.
    !* yges          real       guess radiances.
    !* xn1           real       initial (guess) profile.
    !* np            real       number of retrieved parameters.
    !* nchan         real       number of channels.  
    !*
    !* Outputs:
    !* Var_Name      Type       Description
    !* --------      ----       -----------
    !* xOff          real       offset of linear inversion
    !* xGain         real       gain matrix of linear inversion
    !*
    !* Externals: 
    !* Name    Description
    !* ----    -----------
    !* mxd     multiplication of a matrix with each element of a vector
    !* dminv   matrix inversion
    !*
    !*<copyright>=========================================================
    !*
    !* Developed by Atmospheric and Environmental Research, Inc.                
    !*
    !* Copyright: Atmospheric and Environmental Inc., 1997-2004
    !*                 All Rights Reserved
    !========================================================</copyright>
    !***********************************************************************

    !----Input Variables
    integer,                    intent(in)    :: np, nchan
    real,  dimension(np,nchan), intent(in)    :: xkt
    real,  dimension(nchan),    intent(in)    :: e,yges
    real,  dimension(np,np),    intent(in)    :: cinv
    real,  dimension(np),       intent(in)    :: xn1
    real,  dimension(np),       intent(inout) :: xOff
    real,  dimension(np,nchan), intent(inout) :: xGain
    
    !----Variables of user-defined precision
    real(kind=DBL), dimension(np,nchan)     :: dxkt,dwk,ddy
    real(kind=DBL), dimension(nchan)        :: de,dyges
    real(kind=DBL), dimension(np,np)        :: dcinv,da
    real(kind=DBL), dimension(np)           :: dxn1,dxn,dxOff
    real(kind=DBL), dimension(nchan,np)     :: dxk
    real(kind=DBL), dimension(np,np)        :: dh
    real(kind=DBL)                          :: ddet
    
    !----Conversion from real to user-defined real(kind=DBL)
    dxkt=real(xkt,kind=DBL)
    dcinv=real(cinv,kind=DBL)
    de=real(1./e,kind=DBL)
    dxn1=real(xn1,kind=DBL)
    dyges=real(yges,kind=DBL)

    !----Transpose k matrix
    dxk=transpose(dxkt) 
    !--   calculate contribution functions dy=h-1ktw as (ktwk+g)-1ktw
    call dmxd(dxkt,de,dwk,np,nchan)
    dh = matmul(dwk,dxk)+dcinv
    call dminv(dh,np,ddet)
    ddy=matmul(dh,dwk)
    !----Calculate averaging kernel a=dy*k
    da=matmul(ddy, dxk)
    !---Linear inversion factors
    dxOff=matmul(da, dxn1)-matmul(ddy, dyges)

    !----Conversion from real(kind=DBL) to real
    xOff=real(dxOff)
    xGain=real(ddy)
    
    return
  end subroutine invrt1_GJ_factors


  subroutine invrtRegr(yin,&
       xbak,rbak,wght,uprof,urad,&
       regrcoef,nchan,npg,neof,np,xout)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   invrtRegr
!
! PURPOSE:
!
!   Performs profile inversion by regression.
!
! SYNTAX:
!
!   CALL invrtRegr(yin, xbak, rbak, wght, uprof, urad, regrcoef, 
!      nchan, npg, neof, np, xout)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   yin       REAL     Radiometric data for input
!   xbak      REAL     background state vector
!   rbak      REAL     radiance background
!   wght      REAL     weights
!   uprof     REAL     TBD
!   urad      REAL     TBD
!   regrcoef  REAL     Regression coefficients
!   nchan     INTEGER  Number of channels
!   npg       INTEGER  TBD
!   neof      INTEGER  TBD
!   np        INTEGER  number of parameters
!   
!   INPUTS/OUTPUTS:
!   
!   xout      REAL     TBD
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>


    !----Input Variables
    integer,                      intent(in)    :: np, nchan, neof
    integer,                      intent(in)    :: npg
    real,  dimension(nchan,neof), intent(in)    :: urad
    real,  dimension(npg,np),     intent(in)    :: uprof
    real,  dimension(nchan),      intent(in)    :: yin
    real,  dimension(nchan),      intent(in)    :: wght,rbak
    real,  dimension(neof,np),    intent(in)    :: regrcoef
    real,  dimension(npg),        intent(in)    :: xbak
    real,  dimension(npg),        intent(inout) :: xout
    
    !----Variables of user-defined precision
    real(kind=dbl), dimension(nchan)         :: ytmp
    real(kind=dbl), dimension(neof)          :: ytrans
    real(kind=dbl), dimension(np)            :: dxt
    real(kind=dbl), dimension(npg)           :: dx

    ytmp=yin/wght-rbak
    ytrans=matmul(ytmp,(urad))
    dxt=matmul(ytrans,regrcoef)
    dx=matmul(dxt,transpose(uprof))
    xout=dx+xbak
    return
  end subroutine invrtRegr

  SUBROUTINE dminv(aa,n,d)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   dminv
!
! PURPOSE:
!
!   Inverts a double precision matrix
!
! SYNTAX:
!
!   CALL dminv(aa, n, d)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   n   INTEGER  Dimension of array
!   
!   INPUTS/OUTPUTS:
!   
!   aa  REAL     square matrix
!   d   REAL     Matrix determinant
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

  !--------------------------------------------------------------
  !    * Function name: dminv
  !    * Purpose: Inverts a double precision matrix
  !    * Usage: call dminv(aa,n,d)
  !    * Description: Routine to invert a double precision matrix.
  !    *              The standard gauss-jordan method is used. 
  !    *              the determinant is also calculated. a determinant 
  !    *              of zero indicates that the matrix is singular.
  !    *
  !    * Inputs:
  !    * Var_name   Type       Description
  !    * --------   ----       -----------
  !    * aa         real*8     input matrix (destroyed in computation and 
  !    *                       replaced by resultant inverse.)
  !    * n          integer    order of matrix aa
  !    *
  !    * Outputs:
  !    * Var_name   Type       Description
  !    * --------   ----       -----------
  !    * aa          real*8     resultant inverse matrix
  !    *
  !    * Common blocks: none
  !    * Includes: none
  !    * Externals: none
  !    *
  !    * Copyright: n/a
  !    * Developed by J.R.Eyre
  !    * Modified, March, 2003 by W. Gallery, AER, Inc.
  !    *   Converted to f90 syntax
  !--------------------------------------------------------------

    !---Arguments:
    REAL(KIND=DBL), DIMENSION(n,n), INTENT(INOUT) :: aa
    INTEGER,                        INTENT(IN)    :: n
    REAL(KIND=DBL),                 INTENT(INOUT) :: d

    !---Local variables
    REAL(KIND=DBL), DIMENSION(n)                  :: l,m !previously in arg list
    REAL(KIND=DBL)                                :: biga, hold
    REAL(KIND=DBL), DIMENSION(n*n)                :: a !copy of aa but rank 1

    INTEGER :: i, ij, ik, iz
    INTEGER :: j, ji, jk, jp, jq, jr
    INTEGER :: k, ki, kj, kk
    INTEGER :: nk


    !---Code, as written, reqires rank 1 vector
    a=reshape(aa,(/n*n/))     

    !---search for largest element   
    d=1.0
    nk=-n

    DO k=1,n
       nk=nk+n
       l(k)=k
       m(k)=k
       kk=nk+k
       biga=a(kk)
       DO j=k,n
          iz=n*(j-1)
          DO i=k,n
             ij=iz+i
             IF(abs(biga)-abs(a(ij)) < 0.0) THEN
                biga=a(ij)
                l(k)=i
                m(k)=j
             END IF
          END DO
       END DO
    
       !---interchange rows    
       j=l(k)
       IF((j-k) > 0.0) THEN
          ki=k-n
          DO i=1,n
             ki=ki+n
             hold=-a(ki)
             ji=ki-k+j
             a(ki)=a(ji)
             a(ji) =hold
          END DO
       END IF
   
       !---interchange columns    
       i=m(k)
       IF((i-k) > 0.0) THEN
          jp=n*(i-1)
          DO j=1,n
             jk=nk+j
             ji=jp+j
             hold=-a(jk)
             a(jk)=a(ji)
             a(ji) =hold
          END DO
       END IF
    
       !---divide column by minus pivot (value of pivot element is
       !---contained in biga)    
       IF(biga == 0.0) THEN 
          d=0.0
          RETURN     
       END IF

       DO i=1,n
          IF((i-k) /= 0) THEN
             ik=nk+i
             a(ik)=a(ik)/(-biga)
          END IF
       END DO
             
       !---reduce matrix    
       DO i=1,n
          ik=nk+i
          hold=a(ik)
          ij=i-n
          DO j=1,n
             ij=ij+n
             IF((i-k) /= 0 .and. (j-k) /= 0) THEN
                kj=ij-i+k
                a(ij)=hold*a(kj)+a(ij)
             END IF
          END DO
       END DO    

       !---divide row by pivot    
       kj=k-n
       DO j=1,n
          kj=kj+n
          IF((j-k) /= 0) THEN
             a(kj)=a(kj)/biga
          END IF
       END DO
         
       !---product of pivots   
       !d=d*biga    

       !---replace pivot by reciprocal   
       a(kk)=1.0/biga
    END DO
   
    !---final row and column interchange  
    k=n
    DO 
       k=k-1
       IF(k <= 0) exit
       i=l(k)
       IF((i-k) > 0) THEN 
          jq=n*(k-1)
          jr=n*(i-1)
          DO j=1,n
             jk=jq+j
             hold=a(jk)
             ji=jr+j
             a(jk)=-a(ji)
             a(ji) =hold
          END DO
       END IF
       j=m(k)
       IF((j-k) > 0) THEN
          ki=k-n
          DO i=1,n
             ki=ki+n
             hold=a(ki)
             ji=ki-k+j
             a(ki)=-a(ji)
             a(ji)=hold
          END DO
       END IF
    END DO

    !---Return rank 2 matrix
    aa=reshape(a,(/n,n/))

    RETURN
  END SUBROUTINE dminv
  
  subroutine dmxd(a,b,bcov,nrows,ncols)

!<f90Subroutine>********************************************************
!
! NAME:
!
!   dmxd
!
! PURPOSE:
!
!   Product of matrix by diagonal matrix (stored as a vector). Result is matrix.
!
! SYNTAX:
!
!   CALL dmxd(a, b, bcov, nrows, ncols)
!
! ARGUMENTS:
!
!   INPUTS:
!   
!   a     REAL     Generic array
!   b     REAL     Generic array
!   nrows INTEGER  Number of rows in matrix 
!   ncols INTEGER  Number of columns in matrix 
!   
!   INPUTS/OUTPUTS:
!   
!   bcov  REAL     Background error covariance matrix
!
!   * OPTIONAL
!
! INCLUDES:
!
!   None
!
!*******************************************************</f90Subroutine>

    !*************************************************************************
    !* Function Name:  dmxd
    !* Purpose: Specialized Matrix multiplication
    !* Usage: call dmxd(a,b,bcov,nrows,ncols)
    !* Description: Each column of matrix bcov is the column of matrix a
    !*      multiplied by the corresponding element of row vector b
    !* Inputs:
    !* Var_name     Type       Description
    !* --------     ----       -----------
    !* a            real*8     a nrows by ncols matrix
    !* b            real*8     a vector vector
    !* nrows        integer    row length of the matrix
    !* ncols        integer    column length of the matrix or
    !*                         the length of the vector
    !*
    !* Outputs:
    !* Var_name     Type       Description
    !* --------     ----       -----------
    !* bcov         real *8    product of the matrix and the vector
    !*
    !* Common Blocks: none
    !* Includes: none
    !* Externals: none
    !* Copyright: Atmospheric and Environmental Research, Inc., 1999        
    !* Developed by Atmospheric and Environmental Research, Inc.            
    !    * Modified, March, 2003 by W. Gallery, AER, Inc.
    !    * Converted to f90 syntax
    !***********************************************************************

    integer,                                intent(in)    :: nrows,ncols
    real(kind=DBL), dimension(nrows,ncols), intent(in)    :: a
    real(kind=DBL), dimension(nrows,ncols), intent(inout) :: bcov
    real(kind=DBL), dimension(ncols),       intent(in)    :: b
    integer                          :: i,j

    do j=1,ncols
       do i=1, nrows
          bcov(i,j)=a(i,j)*b(j)
       end do
    end do

  end subroutine dmxd

END MODULE InvertModule
