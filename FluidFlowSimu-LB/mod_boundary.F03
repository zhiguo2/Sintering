!!!--- Module file ---!!!
!!!-Description:
!!!   Subroutines used in boundary condition treatment.
!!!-Details
!!!   Bounce back (half-way)
!!!   Fixed pressure
!!!-Date
!!!   2018-11-27   Eric      V2.0
!!!-Declaration
!!!   Boundary condition subroutine


MODULE bound_mod
USE     mpi_mod, only: pEXCG
USE   const_mod, only: nd,ST,DBL
USE  domain_mod, only: nx,ny,nz
USE lbmodel_mod, only: nq
IMPLICIT NONE


CONTAINS


!*** bounceBack MPI***!
!-Description
!   Perform bounce-back at solid nodes (half-way).
!   Excute exactly after stream process.
!-Arguments
!   f: operated distribution
SUBROUTINE bounceBack(f)
USE  domain_mod, only: i,j,k,pm,pm_solid
USE lbmodel_mod, only: opp,e,iiq
IMPLICIT NONE

REAL(DBL), INTENT(INOUT), DIMENSION(nq,0:(nx+1),0:(ny+1),0:(nz+1))   :: f

REAL(DBL), DIMENSION(nq)     :: f1_temp
INTEGER  :: i_tmp,j_tmp,k_tmp

  DO k=1,nz
  DO j=1,ny
  DO i=1,nx 
    IF (pm(i,j,k).EQ.pm_solid) THEN
    DO iiq=1,nq-1
      i_tmp = i+e(opp(iiq),1)
      j_tmp = j+e(opp(iiq),2)
      k_tmp = k+e(opp(iiq),3)
      f(opp(iiq),i_tmp,j_tmp,k_tmp) = f(iiq,i,j,k)
    END DO
    END IF
  END DO
  END DO
  END DO

END SUBROUTINE bounceBack


!*** bounceBackB MPI***!
!-Description
!   Bounce-back boundary, identified by normal vector index.
!   Perform after stream.
!-Arguments
!   f:  operated distribution
SUBROUTINE bounceBackB(f)
USE     mpi_mod, only: MPI_nbrs,MPI_PROC_NULL
USE  domain_mod, only: i,j,k,pm,pm_solid
USE lbmodel_mod, only: e,iiq
IMPLICIT NONE

REAL(DBL), INTENT(INOUT), DIMENSION(nq,0:(nx+1),0:(ny+1),0:(nz+1))   :: f

INTEGER       :: x_beg,x_end,x_len,y_beg,y_end,y_len,z_beg,z_end,z_len
INTEGER       :: iib,i_tmp,j_tmp,k_tmp
INTEGER(ST), DIMENSION(nq)   :: ops    ! opposite index of velocity

  ops = (/(iiq,iiq=1,nq)/)
  x_len = nx+2
  y_len = ny+2
  z_len = nz+2
  DO iib=1,6       ! boundary index
    x_beg = 0
    x_end = nx+1
    y_beg = 0
    y_end = ny+1
    z_beg = 0
    z_end = nz+1
    SELECT CASE (iib)
    CASE (1)
      x_end = 0
      ops = (/15,15,15, 1,15,15,15,13,15,11,15, 9,15, 7,15/)
    CASE (2)
      y_end = 0
      ops = (/15,15,15,15, 2,15,15,15,12,11,15,15, 8, 7,15/)
    CASE (3)
      z_end = 0
      ops = (/15,15,15,15,15, 3,15,15,15,15,10, 9, 8, 7,15/)
    CASE (4)
      x_beg = nx+1
      ops = (/ 4,15,15,15,15,15,14,15,12,15,10,15, 8,15,15/)
    CASE (5)
      y_beg = ny+1
      ops = (/15, 5,15,15,15,15,14,13,15,15,10, 9,15,15,15/)
    CASE (6)
      z_beg = nz+1
      ops = (/15,15, 6,15,15,15,14,13,12,11,15,15,15,15,15/)
    !CASE DEFAULT
    END SELECT

    DO k=z_beg,z_end
    DO j=y_beg,y_end
    DO i=x_beg,x_end
    IF (pm(i,j,k).EQ.pm_solid) THEN
      DO iiq=1,nq-1
        i_tmp = MOD(x_len+i+e(ops(iiq),1),x_len)
        j_tmp = MOD(y_len+j+e(ops(iiq),2),y_len)
        k_tmp = MOD(z_len+k+e(ops(iiq),3),z_len)
        f(ops(iiq),i_tmp,j_tmp,k_tmp) = f(iiq,i,j,k)
      END DO
    END IF
    END DO
    END DO
    END DO
  END DO

END SUBROUTINE bounceBackB


!*** symmetry MPI ***!
!-Description
!   Symmetry boundary, identified by normal vector index.
!   Perform before stream.
!-Arguments
!   f:  operated distribution
!   bI: boundary index
SUBROUTINE symmetry(f,bI)
USE     mpi_mod, only: MPI_cdnt,MPI_dims
USE  domain_mod, only: i,j,k
USE lbmodel_mod, only: e,iiq
IMPLICIT NONE

REAL(DBL),  INTENT(INOUT), DIMENSION(nq,0:(nx+1),0:(ny+1),0:(nz+1))  :: f
INTEGER,    INTENT(IN)                                               :: bI

INTEGER(ST), DIMENSION(nd)   :: v_nb
INTEGER(ST), DIMENSION(nq)   :: mir    ! mirror index of velocity
INTEGER       :: x_beg,x_end,y_beg,y_end,z_beg,z_end
INTEGER       :: i_tmp,j_tmp,k_tmp
LOGICAL       :: proc_has_bound

  mir = (/(iiq,iiq=1,nq)/)
  proc_has_bound = .FALSE.
  v_nb = e(bI,:)   ! boundary normal vector

  x_beg = 0
  y_beg = 0
  z_beg = 0
  x_end = nx+1
  y_end = ny+1
  z_end = nz+1

  SELECT CASE (bI)
  CASE (1)
    x_end = 0
    mir = (/ 4, 2, 3, 1, 5, 6, 8, 7,10, 9,12,11,14,13,15/)
    IF (MPI_cdnt(1).EQ.0) proc_has_bound = .TRUE.
  CASE (2)
    y_end = 0
    mir = (/ 1, 5, 3, 4, 2, 6, 9,10, 7, 8,13,14,11,12,15/)
    IF (MPI_cdnt(2).EQ.0) proc_has_bound = .TRUE.
  CASE (3)
    z_end = 0
    mir = (/ 1, 2, 6, 4, 5, 3,11,12,13,14, 7, 8, 9,10,15/)
    IF (MPI_cdnt(3).EQ.0) proc_has_bound = .TRUE.
  CASE (4)
    x_beg = nx+1
    mir = (/ 4, 2, 3, 1, 5, 6, 8, 7,10, 9,12,11,14,13,15/)
    IF (MPI_cdnt(1).EQ.(MPI_dims(1)-1)) proc_has_bound = .TRUE.
  CASE (5)
    y_beg = ny+1
    mir = (/ 1, 5, 3, 4, 2, 6, 9,10, 7, 8,13,14,11,12,15/)
    IF (MPI_cdnt(2).EQ.(MPI_dims(2)-1)) proc_has_bound = .TRUE.
  CASE (6)
    z_beg = nz+1
    mir = (/ 1, 2, 6, 4, 5, 3,11,12,13,14, 7, 8, 9,10,15/)
    IF (MPI_cdnt(3).EQ.(MPI_dims(3)-1)) proc_has_bound = .TRUE.
  !CASE DEFAULT
  END SELECT

  IF (proc_has_bound) THEN
  DO k=z_beg,z_end
  DO j=y_beg,y_end
  DO i=x_beg,x_end
  DO iiq=1,nq-1
    i_tmp = i+v_nb(1)
    j_tmp = j+v_nb(2)
    k_tmp = k+v_nb(3)
    f(mir(iiq),i,j,k) = f(iiq,i_tmp,j_tmp,k_tmp)
  END DO
  END DO
  END DO
  END DO
  END IF

END SUBROUTINE symmetry


!*** inletVelocity MPI ***!
!-Description
!   Inlet velocity boundary, identified by normal vector index.
!   Perform beform stream.
!-Arguments
!   um1_in: inlet velocity magnitude
!   bI:     boundary index
SUBROUTINE inletVelocity(um1_in,bI)
USE     mpi_mod, only: MPI_cdnt,MPI_dims
USE  domain_mod, only: i,j,k
USE lbmodel_mod, only: e,computeEqlm
USE    flow_mod, only: f,rho,u
IMPLICIT NONE

REAL(DBL),  INTENT(IN)                                     :: um1_in
INTEGER    ,INTENT(IN)                                     :: bI

INTEGER(ST),DIMENSION(nd)    :: v_nb
REAL(DBL),  DIMENSION(nq)    :: feq1,feq1_n
REAL(DBL),  DIMENSION(nd)    :: u1_in
INTEGER       :: x_beg,x_end,y_beg,y_end,z_beg,z_end
INTEGER       :: i_tmp,j_tmp,k_tmp
LOGICAL       :: proc_has_bound

  proc_has_bound = .FALSE.
  v_nb = e(bI,:)   ! boundary normal vector
  u1_in = um1_in*DBLE(v_nb)

  x_beg = 0
  y_beg = 0
  z_beg = 0
  x_end = nx+1
  y_end = ny+1
  z_end = nz+1

  SELECT CASE (bI)
  CASE (1)
    x_end = 0
    IF (MPI_cdnt(1).EQ.0) proc_has_bound = .TRUE.
  CASE (2)
    y_end = 0
    IF (MPI_cdnt(2).EQ.0) proc_has_bound = .TRUE.
  CASE (3)
    z_end = 0
    IF (MPI_cdnt(3).EQ.0) proc_has_bound = .TRUE.
  CASE (4)
    x_beg = nx+1
    IF (MPI_cdnt(1).EQ.(MPI_dims(1)-1)) proc_has_bound = .TRUE.
  CASE (5)
    y_beg = ny+1
    IF (MPI_cdnt(2).EQ.(MPI_dims(2)-1)) proc_has_bound = .TRUE.
  CASE (6)
    z_beg = nz+1
    IF (MPI_cdnt(3).EQ.(MPI_dims(3)-1)) proc_has_bound = .TRUE.
  !CASE DEFAULT
  END SELECT

  IF (proc_has_bound) THEN
  DO k=z_beg,z_end
  DO j=y_beg,y_end
  DO i=x_beg,x_end
    i_tmp = i+v_nb(1)
    j_tmp = j+v_nb(2)
    k_tmp = k+v_nb(3)
    CALL computeEqlm(feq1  ,rho(i_tmp,j_tmp,k_tmp),u1_in                 )
    CALL computeEqlm(feq1_n,rho(i_tmp,j_tmp,k_tmp),u(:,i_tmp,j_tmp,k_tmp))
    f(:,i,j,k) = feq1 + (f(:,i_tmp,j_tmp,k_tmp) - feq1_n)
  END DO
  END DO
  END DO
  END IF

END SUBROUTINE inletVelocity


!*** fixedPressure MPI ***!
!-Description
!   Fixed pressure boundary, identified by normal vector index.
!   Perform beform stream.
!-Arguments
!   dlt_p1: given pressure difference
!   bI:     boundary index
SUBROUTINE fixedPressure(dlt_p1,bI)
USE     mpi_mod, only: MPI_cdnt,MPI_dims
USE  domain_mod, only: i,j,k
USE lbmodel_mod, only: e,RT,computeEqlm
USE    flow_mod, only: f,rho,u,rho0
IMPLICIT NONE

REAL(DBL), INTENT(IN)                                      :: dlt_p1
INTEGER,   INTENT(IN)                                      :: bI

INTEGER(ST),DIMENSION(nd)    :: v_nb
REAL(DBL),  DIMENSION(nq)    :: feq1,feq1_n
INTEGER       :: x_beg,x_end,y_beg,y_end,z_beg,z_end
INTEGER       :: i_tmp,j_tmp,k_tmp
REAL(DBL)     :: rho_temp
LOGICAL       :: proc_has_bound

  proc_has_bound = .FALSE.
  v_nb = e(bI,:)   ! boundary normal vector

  x_beg = 0
  y_beg = 0
  z_beg = 0
  x_end = nx+1
  y_end = ny+1
  z_end = nz+1

  SELECT CASE (bI)
  CASE (1)
    x_end = 0
    IF (MPI_cdnt(1).EQ.0) proc_has_bound = .TRUE.
  CASE (2)
    y_end = 0
    IF (MPI_cdnt(2).EQ.0) proc_has_bound = .TRUE.
  CASE (3)
    z_end = 0
    IF (MPI_cdnt(3).EQ.0) proc_has_bound = .TRUE.
  CASE (4)
    x_beg = nx+1
    IF (MPI_cdnt(1).EQ.(MPI_dims(1)-1)) proc_has_bound = .TRUE.
  CASE (5)
    y_beg = ny+1
    IF (MPI_cdnt(2).EQ.(MPI_dims(2)-1)) proc_has_bound = .TRUE.
  CASE (6)
    z_beg = nz+1
    IF (MPI_cdnt(3).EQ.(MPI_dims(3)-1)) proc_has_bound = .TRUE.
  !CASE DEFAULT
  END SELECT

  IF (proc_has_bound) THEN
  DO k=z_beg,z_end
  DO j=y_beg,y_end
  DO i=x_beg,x_end
    i_tmp = i+v_nb(1)
    j_tmp = j+v_nb(2)
    k_tmp = k+v_nb(3)
    !IF (DOT_PRODUCT(u(:,i_tmp,j_tmp,k_tmp),v_nb).LT.0) THEN
    !  rho_temp = rho(i_tmp,j_tmp,k_tmp)
    !ELSE
      rho_temp = rho0 + dlt_p1/RT
    !END IF
    CALL computeEqlm(feq1  ,rho_temp              ,u(:,i_tmp,j_tmp,k_tmp))
    CALL computeEqlm(feq1_n,rho(i_tmp,j_tmp,k_tmp),u(:,i_tmp,j_tmp,k_tmp))
    f(:,i,j,k) = feq1 + (f(:,i_tmp,j_tmp,k_tmp) - feq1_n)*rho_temp/rho(i_tmp,j_tmp,k_tmp)
  END DO
  END DO
  END DO
  END IF

END SUBROUTINE fixedPressure

END MODULE bound_mod
