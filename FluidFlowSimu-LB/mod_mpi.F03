!!!--- Module file ---!!!
!!!-Description:
!!!    MPI related information.
!!!-Details
!!!   MPI parameters
!!!   MPI virtual topology related
!!!   MPI reduce subroutines
!!!-Date
!!!   2018-11-27   Eric      V2.0
!!!-Declaration
!!!   MPI module


MODULE mpi_mod
USE   mpi
USE const_mod, only: nd,ST,DBL
IMPLICIT NONE

!*** Generic reduce functions ***! 
! Reduce maxval
INTERFACE pMAXVAL
  MODULE PROCEDURE MAXREAL3,MAXREAL4
END INTERFACE pMAXVAL

! Reduce minval
INTERFACE pMINVAL
  MODULE PROCEDURE MINREAL3,MINREAL4
END INTERFACE pMINVAL

! Reduce sum
INTERFACE pSUM
  MODULE PROCEDURE SUMINT3,SUMINT4,SUMREAL3,SUMREAL4
END INTERFACE pSUM


!*** Generic message functions ***!
INTERFACE pEXCG
  MODULE PROCEDURE EXCGSINT3,EXCGREAL3,EXCGREAL4,EXCGREAL4E
END INTERFACE pEXCG


!*** MPI parameters ***!
INTEGER  :: MPI_ierr         ! error identifier
INTEGER  :: MPI_id,MPI_np    ! current proc id, number of procs
INTEGER  :: iip              ! process iter index
INTEGER  :: MPI_comm_new     ! new communication
INTEGER  :: MPI_stat(MPI_STATUS_SIZE)  ! status


!*** MPI Vitual Topology parameters ***!
LOGICAL  :: MPI_reod         ! reorder in new comm
LOGICAL, DIMENSION(nd)   :: MPI_perd
INTEGER, DIMENSION(nd)   :: MPI_dims, MPI_cdnt
                        ! dimension size, whether periodic, coordinate in new comm
INTEGER, DIMENSION(nd*2) :: MPI_nbrs ! neighbours - E N W S
INTEGER                  :: MPI_dims_idct ! dims indictor 0 - create dims automatically


CONTAINS


! Obtain subdir
!-Description
!   Get output dir for proc
!-Arguments
!   iip: proc id
FUNCTION subdir(iip) RESULT(res)
IMPLICIT NONE

INTEGER, INTENT(IN)     :: iip
CHARACTER(11)           :: res

CHARACTER(2)  :: procI

  WRITE(procI,'(I2.2)') iip
  res = 'processor'//TRIM(procI)

END FUNCTION subdir


! Initialize MPI processes and parameters
SUBROUTINE initMPI
IMPLICIT NONE

CHARACTER(70) :: dummy
INTEGER(ST), DIMENSION(2*nd) :: bound

  ! Input parameter
  OPEN (UNIT=41,FILE='INPUT/mpiDict',STATUS='UNKNOWN')
  READ(41,*)
  READ(41,*) dummy,MPI_dims_idct
  IF (MPI_dims_idct.NE.0) THEN
    READ(41,*) dummy,MPI_dims
  END IF
  CLOSE (UNIT=41)

  ! Whether periodic - topology relation
  OPEN (UNIT=41,FILE='INPUT/domainDict',STATUS='UNKNOWN')
  READ(41,*)
  READ(41,*)
  READ(41,*)
  READ(41,*)
  READ(41,*)
  READ(41,*)
  READ(41,*) dummy,bound
  CLOSE (UNIT=41)

  MPI_reod = .TRUE.
  MPI_perd = (/.FALSE.,.FALSE.,.FALSE./)
  IF (bound(1).EQ.0 .AND. bound(4).EQ.0) MPI_perd(1) = .TRUE.
  IF (bound(2).EQ.0 .AND. bound(5).EQ.0) MPI_perd(2) = .TRUE.
  IF (bound(3).EQ.0 .AND. bound(6).EQ.0) MPI_perd(3) = .TRUE.

  ! Initialize MPI process
  CALL MPI_INIT(MPI_ierr)
  ! Get proc id, number of procs
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MPI_id, MPI_ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, MPI_np, MPI_ierr)
  ! Create Cartesian topology
  IF (MPI_dims_idct.EQ.0) THEN
  CALL MPI_DIMS_CREATE(MPI_np, nd, MPI_dims, MPI_ierr)
  END IF
  CALL MPI_CART_CREATE(MPI_COMM_WORLD, nd, MPI_dims, MPI_perd, MPI_reod, MPI_comm_new, MPI_ierr)
  ! Get coordinate and neighbours
  CALL MPI_CART_COORDS(MPI_comm_new, MPI_id, nd, MPI_cdnt, MPI_ierr)
  CALL MPI_CART_SHIFT( MPI_comm_new, 0, 1, MPI_nbrs(4), MPI_nbrs(1), MPI_ierr)
  CALL MPI_CART_SHIFT( MPI_comm_new, 1, 1, MPI_nbrs(5), MPI_nbrs(2), MPI_ierr)
  CALL MPI_CART_SHIFT( MPI_comm_new, 2, 1, MPI_nbrs(6), MPI_nbrs(3), MPI_ierr)

  ! Create output directory if needed
  CALL SYSTEM('mkdir -p OUTPUT/'//TRIM(subdir(MPI_id)))

  ! Print
  CALL printMPI

END SUBROUTINE initMPI


! print MPI processes relations
SUBROUTINE printMPI
IMPLICIT NONE

  ! write MPI information
  10 FORMAT (A)
  11 FORMAT (I3,1X,A,6I5)
  IF (MPI_id.EQ.0) THEN
    WRITE(*,10) '********* Create MPI processes ********'
    WRITE(*,10) ' ID |    R    F    U    L    B    D    '
    WRITE(*,10) '---------------------------------------'
  END IF
  DO iip=0,MPI_np-1
    CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_ierr)
    IF (MPI_id.EQ.iip) THEN
      WRITE(*,11) MPI_id,'|',MPI_nbrs
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_ierr)
  END DO
  CALL MPI_BARRIER(MPI_COMM_WORLD,MPI_ierr)
  IF (MPI_id.EQ.0) THEN
    WRITE(*,10) '---------------------------------------'
    WRITE(*,10) '** Successfully create MPI processes **'
    WRITE(*,*)
  END IF

END SUBROUTINE printMPI


! Specific reduce functions
! pMAXVAL
FUNCTION MAXREAL3(arr_real) RESULT(max_val)
IMPLICIT NONE

REAL(DBL), INTENT(IN), DIMENSION(:,:,:)     :: arr_real
REAL(DBL)     :: mymax,max_val

  mymax = MAXVAL(arr_real)
  CALL MPI_REDUCE(mymax,max_val,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPI_ierr)

END FUNCTION MAXREAL3

FUNCTION MAXREAL4(arr_real) RESULT(max_val)
IMPLICIT NONE

REAL(DBL), INTENT(IN), DIMENSION(:,:,:,:)   :: arr_real
REAL(DBL)     :: mymax,max_val

  mymax = MAXVAL(arr_real)
  CALL MPI_REDUCE(mymax,max_val,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,MPI_ierr)

END FUNCTION MAXREAL4


! pMINVAL
FUNCTION MINREAL3(arr_real) RESULT(min_val)
IMPLICIT NONE

REAL(DBL), INTENT(IN), DIMENSION(:,:,:)     :: arr_real
REAL(DBL)     :: mymin,min_val

  mymin = MINVAL(arr_real)
  CALL MPI_REDUCE(mymin,min_val,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,MPI_ierr)

END FUNCTION MINREAL3

FUNCTION MINREAL4(arr_real) RESULT(min_val)
IMPLICIT NONE

REAL(DBL), INTENT(IN), DIMENSION(:,:,:,:)   :: arr_real
REAL(DBL)     :: mymin,min_val

  mymin = MINVAL(arr_real)
  CALL MPI_REDUCE(mymin,min_val,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,MPI_ierr)

END FUNCTION MINREAL4


! pSUM
FUNCTION SUMINT3(arr_int) RESULT(sum_int)
IMPLICIT NONE
INTEGER, INTENT(IN), DIMENSION(:,:,:)       :: arr_int
INTEGER       :: mysum,sum_int

  mysum = SUM(arr_int)
  CALL MPI_REDUCE(mysum,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,MPI_ierr)

END FUNCTION SUMINT3

FUNCTION SUMINT4(arr_int) RESULT(sum_int)
IMPLICIT NONE
INTEGER, INTENT(IN), DIMENSION(:,:,:,:)     :: arr_int
INTEGER       :: mysum,sum_int

  mysum = SUM(arr_int)
  CALL MPI_REDUCE(mysum,sum_int,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,MPI_ierr)

END FUNCTION SUMINT4

FUNCTION SUMREAL3(arr_real) RESULT(sum_real)
IMPLICIT NONE

REAL(DBL), INTENT(IN), DIMENSION(:,:,:)     :: arr_real
REAL(DBL)     :: mysum,sum_real

  mysum = SUM(arr_real)
  CALL MPI_REDUCE(mysum,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPI_ierr)

END FUNCTION SUMREAL3

FUNCTION SUMREAL4(arr_real) RESULT(sum_real)
IMPLICIT NONE

REAL(DBL), INTENT(IN), DIMENSION(:,:,:,:)   :: arr_real
REAL(DBL)     :: mysum,sum_real

  mysum = SUM(arr_real)
  CALL MPI_REDUCE(mysum,sum_real,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPI_ierr)

END FUNCTION SUMREAL4


!pEXCG
SUBROUTINE EXCGSINT3(a)
IMPLICIT NONE

INTEGER(ST), INTENT(INOUT), DIMENSION(0:,0:,0:)  :: a

INTEGER,     ALLOCATABLE,   DIMENSION(:,:)       :: ax_in,ax_out
INTEGER,     ALLOCATABLE,   DIMENSION(:,:)       :: ay_in,ay_out
INTEGER,     ALLOCATABLE,   DIMENSION(:,:)       :: az_in,az_out
INTEGER       :: nx,ny,nz
INTEGER       :: a_len,next,prev


  nx = SIZE(a,1)-2
  ny = SIZE(a,2)-2
  nz = SIZE(a,3)-2

  IF (.NOT.ALLOCATED(ax_in )) ALLOCATE(ax_in( ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ax_out)) ALLOCATE(ax_out(ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_in )) ALLOCATE(ay_in( nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_out)) ALLOCATE(ay_out(nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_in )) ALLOCATE(az_in( nx+2,ny+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_out)) ALLOCATE(az_out(nx+2,ny+2),STAT=MPI_ierr)

  ! R F U L B D
  ! R - Right
  ax_out = INT(a(nx,:,:))
  ax_in  = INT(a(0 ,:,:))
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(1)
  prev  = MPI_nbrs(4)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_INTEGER,next,1,&
                   &ax_in ,a_len,MPI_INTEGER,prev,1,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(0,:,:) = INT(ax_in,ST)

  ! F - Front
  ay_out = INT(a(:,ny,:))
  ay_in  = INT(a(:,0 ,:))
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(2)
  prev  = MPI_nbrs(5)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_INTEGER,next,2,&
                   &ay_in ,a_len,MPI_INTEGER,prev,2,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,0,:) = INT(ay_in,ST)

  ! U - Up
  az_out = INT(a(:,:,nz))
  az_in  = INT(a(:,:,0 ))
  a_len = SIZE(az_out)
  next  = MPI_nbrs(3)
  prev  = MPI_nbrs(6)
  CALL MPI_SENDRECV(az_out,a_len,MPI_INTEGER,next,3,&
                   &az_in ,a_len,MPI_INTEGER,prev,3,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,0) = INT(az_in,ST)

  ! L - Left
  ax_out = INT(a(   1,:,:))
  ax_in  = INT(a(nx+1,:,:))
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(4)
  prev  = MPI_nbrs(1)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_INTEGER,next,4,&
                   &ax_in ,a_len,MPI_INTEGER,prev,4,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(nx+1,:,:) = INT(ax_in,ST)

  ! B - Back
  ay_out = INT(a(:,   1,:))
  ay_in  = INT(a(:,ny+1,:))
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(5)
  prev  = MPI_nbrs(2)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_INTEGER,next,5,&
                   &ay_in ,a_len,MPI_INTEGER,prev,5,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,ny+1,:) = INT(ay_in,ST)

  ! D - Down
  az_out = INT(a(:,:,   1))
  az_in  = INT(a(:,:,nz+1))
  a_len = SIZE(az_out)
  next  = MPI_nbrs(6)
  prev  = MPI_nbrs(3)
  CALL MPI_SENDRECV(az_out,a_len,MPI_INTEGER,next,6,&
                   &az_in ,a_len,MPI_INTEGER,prev,6,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,nz+1) = INT(az_in)

END SUBROUTINE EXCGSINT3

SUBROUTINE EXCGREAL3(a)
IMPLICIT NONE

REAL(DBL), INTENT(INOUT), DIMENSION(0:,0:,0:)    :: a

REAL(DBL), ALLOCATABLE,   DIMENSION(:,:)         :: ax_in,ax_out
REAL(DBL), ALLOCATABLE,   DIMENSION(:,:)         :: ay_in,ay_out
REAL(DBL), ALLOCATABLE,   DIMENSION(:,:)         :: az_in,az_out
INTEGER       :: nx,ny,nz
INTEGER       :: a_len,next,prev


  nx = SIZE(a,1)-2
  ny = SIZE(a,2)-2
  nz = SIZE(a,3)-2

  IF (.NOT.ALLOCATED(ax_in )) ALLOCATE(ax_in( ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ax_out)) ALLOCATE(ax_out(ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_in )) ALLOCATE(ay_in( nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_out)) ALLOCATE(ay_out(nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_in )) ALLOCATE(az_in( nx+2,ny+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_out)) ALLOCATE(az_out(nx+2,ny+2),STAT=MPI_ierr)

  ! R F U L B D
  ! R - Right
  ax_out = a(nx,:,:)
  ax_in  = a(0 ,:,:)
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(1)
  prev  = MPI_nbrs(4)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_DOUBLE_PRECISION,next,1,&
                   &ax_in ,a_len,MPI_DOUBLE_PRECISION,prev,1,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(0,:,:) = ax_in

  ! F - Front
  ay_out = a(:,ny,:)
  ay_in  = a(:,0 ,:)
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(2)
  prev  = MPI_nbrs(5)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_DOUBLE_PRECISION,next,2,&
                   &ay_in ,a_len,MPI_DOUBLE_PRECISION,prev,2,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,0,:) = ay_in

  ! U - Up
  az_out = a(:,:,nz)
  az_in  = a(:,:,0 )
  a_len = SIZE(az_out)
  next  = MPI_nbrs(3)
  prev  = MPI_nbrs(6)
  CALL MPI_SENDRECV(az_out,a_len,MPI_DOUBLE_PRECISION,next,3,&
                   &az_in ,a_len,MPI_DOUBLE_PRECISION,prev,3,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,0) = az_in

  ! L - Left
  ax_out = a(   1,:,:)
  ax_in  = a(nx+1,:,:)
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(4)
  prev  = MPI_nbrs(1)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_DOUBLE_PRECISION,next,4,&
                   &ax_in ,a_len,MPI_DOUBLE_PRECISION,prev,4,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(nx+1,:,:) = ax_in

  ! B - Back
  ay_out = a(:,   1,:)
  ay_in  = a(:,ny+1,:)
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(5)
  prev  = MPI_nbrs(2)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_DOUBLE_PRECISION,next,5,&
                   &ay_in ,a_len,MPI_DOUBLE_PRECISION,prev,5,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,ny+1,:) = ay_in

  ! D - Down
  az_out = a(:,:,   1)
  az_in  = a(:,:,nz+1)
  a_len = SIZE(az_out)
  next  = MPI_nbrs(6)
  prev  = MPI_nbrs(3)
  CALL MPI_SENDRECV(az_out,a_len,MPI_DOUBLE_PRECISION,next,6,&
                   &az_in ,a_len,MPI_DOUBLE_PRECISION,prev,6,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,nz+1) = az_in

END SUBROUTINE EXCGREAL3

SUBROUTINE EXCGREAL4(a)
IMPLICIT NONE

REAL(DBL), INTENT(INOUT), DIMENSION(:,0:,0:,0:)  :: a

REAL(DBL), ALLOCATABLE,   DIMENSION(:,:,:)       :: ax_in,ax_out
REAL(DBL), ALLOCATABLE,   DIMENSION(:,:,:)       :: ay_in,ay_out
REAL(DBL), ALLOCATABLE,   DIMENSION(:,:,:)       :: az_in,az_out
INTEGER       :: nq,nx,ny,nz
INTEGER       :: a_len,next,prev

  nq = SIZE(a,1)
  nx = SIZE(a,2)-2
  ny = SIZE(a,3)-2
  nz = SIZE(a,4)-2

  IF (.NOT.ALLOCATED(ax_in )) ALLOCATE(ax_in( nq,ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ax_out)) ALLOCATE(ax_out(nq,ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_in )) ALLOCATE(ay_in( nq,nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_out)) ALLOCATE(ay_out(nq,nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_in )) ALLOCATE(az_in( nq,nx+2,ny+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_out)) ALLOCATE(az_out(nq,nx+2,ny+2),STAT=MPI_ierr)

  ! R F U L B D
  ! R - Right
  ax_out = a(:,nx,:,:)
  ax_in  = a(:,0 ,:,:)
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(1)
  prev  = MPI_nbrs(4)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_DOUBLE_PRECISION,next,1,&
                   &ax_in ,a_len,MPI_DOUBLE_PRECISION,prev,1,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,0,:,:) = ax_in

  ! F - Front
  ay_out = a(:,:,ny,:)
  ay_in  = a(:,:,0 ,:)
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(2)
  prev  = MPI_nbrs(5)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_DOUBLE_PRECISION,next,2,&
                   &ay_in ,a_len,MPI_DOUBLE_PRECISION,prev,2,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,0,:) = ay_in

  ! U - Up
  az_out = a(:,:,:,nz)
  az_in  = a(:,:,:,0 )
  a_len = SIZE(az_out)
  next  = MPI_nbrs(3)
  prev  = MPI_nbrs(6)
  CALL MPI_SENDRECV(az_out,a_len,MPI_DOUBLE_PRECISION,next,3,&
                   &az_in ,a_len,MPI_DOUBLE_PRECISION,prev,3,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,:,0) = az_in

  ! L - Left
  ax_out = a(:,   1,:,:)
  ax_in  = a(:,nx+1,:,:)
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(4)
  prev  = MPI_nbrs(1)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_DOUBLE_PRECISION,next,4,&
                   &ax_in ,a_len,MPI_DOUBLE_PRECISION,prev,4,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,nx+1,:,:) = ax_in

  ! B - Back
  ay_out = a(:,:,   1,:)
  ay_in  = a(:,:,ny+1,:)
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(5)
  prev  = MPI_nbrs(2)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_DOUBLE_PRECISION,next,5,&
                   &ay_in ,a_len,MPI_DOUBLE_PRECISION,prev,5,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,ny+1,:) = ay_in

  ! D - Down
  az_out = a(:,:,:,   1)
  az_in  = a(:,:,:,nz+1)
  a_len = SIZE(az_out)
  next  = MPI_nbrs(6)
  prev  = MPI_nbrs(3)
  CALL MPI_SENDRECV(az_out,a_len,MPI_DOUBLE_PRECISION,next,6,&
                   &az_in ,a_len,MPI_DOUBLE_PRECISION,prev,6,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,:,nz+1) = az_in

END SUBROUTINE EXCGREAL4

SUBROUTINE EXCGREAL4E(a,e)
USE const_mod, only: ST,DBL
IMPLICIT NONE

REAL(DBL), INTENT(INOUT), DIMENSION(:,0:,0:,0:)  :: a
INTEGER(ST),INTENT(IN),   DIMENSION(:,:)         :: e

REAL(DBL), ALLOCATABLE,   DIMENSION(:,:)         :: ax_in,ax_out
REAL(DBL), ALLOCATABLE,   DIMENSION(:,:)         :: ay_in,ay_out
REAL(DBL), ALLOCATABLE,   DIMENSION(:,:)         :: az_in,az_out
INTEGER       :: iiq,nq,nx,ny,nz
INTEGER       :: a_len,next,prev

  nq = SIZE(a,1)
  nx = SIZE(a,2)-2
  ny = SIZE(a,3)-2
  nz = SIZE(a,4)-2

  IF (SIZE(e,1).EQ.nq .AND. SIZE(e,2).EQ.3) THEN

  IF (.NOT.ALLOCATED(ax_in )) ALLOCATE(ax_in( ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ax_out)) ALLOCATE(ax_out(ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_in )) ALLOCATE(ay_in( nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_out)) ALLOCATE(ay_out(nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_in )) ALLOCATE(az_in( nx+2,ny+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_out)) ALLOCATE(az_out(nx+2,ny+2),STAT=MPI_ierr)

  ! R F U L B D
  ! R - Right
  DO iiq=1,nq
  IF (e(iiq,1).EQ.1) THEN
  ax_out = a(iiq,nx+1,:,:)
  ax_in  = a(iiq,1 ,:,:)
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(1)
  prev  = MPI_nbrs(4)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_DOUBLE_PRECISION,next,1,&
                   &ax_in ,a_len,MPI_DOUBLE_PRECISION,prev,1,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(iiq,1,:,:) = ax_in
  END IF
  END DO

  ! F - Front
  DO iiq=1,nq
  IF (e(iiq,2).EQ.1) THEN
  ay_out = a(iiq,:,ny+1,:)
  ay_in  = a(iiq,:,1 ,:)
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(2)
  prev  = MPI_nbrs(5)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_DOUBLE_PRECISION,next,2,&
                   &ay_in ,a_len,MPI_DOUBLE_PRECISION,prev,2,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(iiq,:,1,:) = ay_in
  END IF
  END DO

  ! U - Up
  DO iiq=1,nq
  IF (e(iiq,3).EQ.1) THEN
  az_out = a(iiq,:,:,nz+1)
  az_in  = a(iiq,:,:,1 )
  a_len = SIZE(az_out)
  next  = MPI_nbrs(3)
  prev  = MPI_nbrs(6)
  CALL MPI_SENDRECV(az_out,a_len,MPI_DOUBLE_PRECISION,next,3,&
                   &az_in ,a_len,MPI_DOUBLE_PRECISION,prev,3,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(iiq,:,:,1) = az_in
  END IF
  END DO

  ! L - Left
  DO iiq=1,nq
  IF (e(iiq,1).EQ.-1) THEN
  ax_out = a(iiq, 0,:,:)
  ax_in  = a(iiq,nx,:,:)
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(4)
  prev  = MPI_nbrs(1)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_DOUBLE_PRECISION,next,4,&
                   &ax_in ,a_len,MPI_DOUBLE_PRECISION,prev,4,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(iiq,nx,:,:) = ax_in
  END IF
  END DO

  ! B - Back
  DO iiq=1,nq
  IF (e(iiq,2).EQ.-1) THEN
  ay_out = a(iiq,:, 0,:)
  ay_in  = a(iiq,:,ny,:)
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(5)
  prev  = MPI_nbrs(2)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_DOUBLE_PRECISION,next,5,&
                   &ay_in ,a_len,MPI_DOUBLE_PRECISION,prev,5,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(iiq,:,ny,:) = ay_in
  END IF
  END DO

  ! D - Down
  DO iiq=1,nq
  IF (e(iiq,3).NE.1) THEN
  az_out = a(iiq,:,:, 0)
  az_in  = a(iiq,:,:,nz)
  a_len = SIZE(az_out)
  next  = MPI_nbrs(6)
  prev  = MPI_nbrs(3)
  CALL MPI_SENDRECV(az_out,a_len,MPI_DOUBLE_PRECISION,next,6,&
                   &az_in ,a_len,MPI_DOUBLE_PRECISION,prev,6,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(iiq,:,:,nz) = az_in
  END IF
  END DO

  END IF

END SUBROUTINE EXCGREAL4E


! exchange and add
SUBROUTINE pEXCGA(a)
IMPLICIT NONE

REAL(DBL), INTENT(INOUT), DIMENSION(:,0:,0:,0:)  :: a

REAL(DBL), ALLOCATABLE,   DIMENSION(:,:,:)       :: ax_in,ax_out
REAL(DBL), ALLOCATABLE,   DIMENSION(:,:,:)       :: ay_in,ay_out
REAL(DBL), ALLOCATABLE,   DIMENSION(:,:,:)       :: az_in,az_out
INTEGER       :: nq,nx,ny,nz
INTEGER       :: a_len,next,prev

  nq = SIZE(a,1)
  nx = SIZE(a,2)-2
  ny = SIZE(a,3)-2
  nz = SIZE(a,4)-2

  IF (.NOT.ALLOCATED(ax_in )) ALLOCATE(ax_in( nq,ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ax_out)) ALLOCATE(ax_out(nq,ny+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_in )) ALLOCATE(ay_in( nq,nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(ay_out)) ALLOCATE(ay_out(nq,nx+2,nz+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_in )) ALLOCATE(az_in( nq,nx+2,ny+2),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(az_out)) ALLOCATE(az_out(nq,nx+2,ny+2),STAT=MPI_ierr)

  ! R F U L B D
  ! R - Right
  ax_out = a(:,nx+1,:,:)
  ax_in  = 0d0*a(:,   1,:,:)
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(1)
  prev  = MPI_nbrs(4)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_DOUBLE_PRECISION,next,1,&
                   &ax_in ,a_len,MPI_DOUBLE_PRECISION,prev,1,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,1,:,:) = a(:,1,:,:) + ax_in

  ! F - Front
  ay_out = a(:,:,ny+1,:)
  ay_in  = 0d0*a(:,:,   1,:)
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(2)
  prev  = MPI_nbrs(5)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_DOUBLE_PRECISION,next,2,&
                   &ay_in ,a_len,MPI_DOUBLE_PRECISION,prev,2,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,1,:) = a(:,:,1,:) + ay_in

  ! U - Up
  az_out = a(:,:,:,nz+1)
  az_in  = 0d0*a(:,:,:,   1)
  a_len = SIZE(az_out)
  next  = MPI_nbrs(3)
  prev  = MPI_nbrs(6)
  CALL MPI_SENDRECV(az_out,a_len,MPI_DOUBLE_PRECISION,next,3,&
                   &az_in ,a_len,MPI_DOUBLE_PRECISION,prev,3,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,:,1) = a(:,:,:,1) + az_in

  ! L - Left
  ax_out = a(:, 0,:,:)
  ax_in  = 0d0*a(:,nx,:,:)
  a_len = SIZE(ax_out)
  next  = MPI_nbrs(4)
  prev  = MPI_nbrs(1)
  CALL MPI_SENDRECV(ax_out,a_len,MPI_DOUBLE_PRECISION,next,4,&
                   &ax_in ,a_len,MPI_DOUBLE_PRECISION,prev,4,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,nx,:,:) = a(:,nx,:,:) + ax_in

  ! B - Back
  ay_out = a(:,:, 0,:)
  ay_in  = 0d0*a(:,:,ny,:)
  a_len = SIZE(ay_out)
  next  = MPI_nbrs(5)
  prev  = MPI_nbrs(2)
  CALL MPI_SENDRECV(ay_out,a_len,MPI_DOUBLE_PRECISION,next,5,&
                   &ay_in ,a_len,MPI_DOUBLE_PRECISION,prev,5,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,ny,:) = a(:,:,ny,:) + ay_in

  ! D - Down
  az_out = a(:,:,:, 0)
  az_in  = 0d0*a(:,:,:,nz)
  a_len = SIZE(az_out)
  next  = MPI_nbrs(6)
  prev  = MPI_nbrs(3)
  CALL MPI_SENDRECV(az_out,a_len,MPI_DOUBLE_PRECISION,next,6,&
                   &az_in ,a_len,MPI_DOUBLE_PRECISION,prev,6,&
                   &MPI_COMM_WORLD,MPI_stat,MPI_ierr)
  a(:,:,:,nz) = a(:,:,:,nz) + az_in

END SUBROUTINE pEXCGA

END MODULE mpi_mod
