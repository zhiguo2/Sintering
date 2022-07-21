!!!--- Module file ---!!!
!!!-Description:
!!!   Doamin information in both space and time.
!!!-Details
!!!   Domain dimension and size
!!!   Time steps and control
!!!   Solids obstacles
!!!-Date
!!!   2018-11-27   Eric      V2.0
!!!-Declaration
!!!   Domain module


MODULE domain_mod
USE const_mod, only: nd,ST,LG,DBL
USE   mpi_mod
IMPLICIT NONE

!*** Domain size ***!
INTEGER(LG)             :: nx_t,ny_t,nz_t   ! domain size (total)
INTEGER(LG)             :: nx,ny,nz         ! domain size (local proc)
REAL(DBL)               :: lx_t,ly_t,lz_t   ! physical length (total)
REAL(DBL)               :: dx,dy,dz         ! lattice step
REAL(DBL)               :: uc               ! lattice speed
INTEGER(LG)             :: i,j,k            ! space index
INTEGER, DIMENSION(nd)  :: MPI_orig         ! origin coordination
INTEGER, DIMENSION(:,:), ALLOCATABLE        :: MPI_block
                                            ! gather each processor's block extent

!*** Time control ***!
INTEGER(LG)             :: nt               ! total time
INTEGER(LG)             :: nt_s             ! start time
INTEGER(LG)             :: nt_o             ! time interval to output macro properties
INTEGER(LG)             :: nt_f             ! time interval to output distribution functions
REAL(DBL)               :: lt               ! physical time length
REAL(DBL)               :: dt               ! physical time step
INTEGER(LG)             :: iit              ! time index


!*** Geometry ***!
INTEGER(ST), DIMENSION(:,:,:), ALLOCATABLE  :: pm     ! porous media
REAL(DBL), DIMENSION(:,:,:,:), ALLOCATABLE  :: v_nw   ! wall normal vector
INTEGER(ST), DIMENSION(2*nd)                :: bound 
                                            ! boundary index - 0 for cyclic; 1 for bounce back
                                            ! others for other boundary
INTEGER                                     :: iib
CHARACTER(120)          :: pm_file          ! file that stores porous media structure info
INTEGER(ST)             :: pm_solid         ! solid point identifier (usual value is 1)
INTEGER(ST)             :: pm_pore          !  pore point identifier (usual value is 0)

CONTAINS

SUBROUTINE initDomain
IMPLICIT NONE

CHARACTER(70) :: dummy
INTEGER       :: int_tmp,nc  ! characteristic 
REAL(DBL)     :: lc

  ! Input
  OPEN (UNIT=41,FILE='INPUT/domainDict',STATUS='UNKNOWN')
  READ(41,*)
  READ(41,*) dummy,nx_t,ny_t,nz_t
  READ(41,*) dummy,nt,nt_o,nt_f,nt_s
  READ(41,*) dummy,nc,lc
  READ(41,*) dummy,pm_file
  READ(41,*) dummy,pm_pore,pm_solid
  READ(41,*) dummy,bound
  CLOSE (UNIT=41)
  dx = lc / nc; dy = dx; dz = dx
  lx_t = nx_t * dx
  ly_t = ny_t * dy
  lz_t = nz_t * dz

  ! Get local size 
  CALL setSize

  ! Allocate memory
  IF (.NOT.ALLOCATED(  pm)) ALLOCATE(  pm(   0:(nx+1),0:(ny+1),0:(nz+1)),STAT=MPI_ierr)

  ! Initialize solid
  pm = pm_pore
  CALL readPm

  ! Gather block extent
  IF (MPI_id.EQ.0) THEN
    IF (.NOT.ALLOCATED(MPI_block)) ALLOCATE(MPI_block(MPI_np,2*nd),STAT=MPI_ierr)
    MPI_block = 0
  END IF
  int_tmp = MPI_orig(1)
  CALL MPI_GATHER(int_tmp,1,MPI_INT,MPI_block(1,1),1,MPI_INT,0,MPI_COMM_WORLD,MPI_ierr)
  int_tmp = MPI_orig(1)+nx
  CALL MPI_GATHER(int_tmp,1,MPI_INT,MPI_block(1,2),1,MPI_INT,0,MPI_COMM_WORLD,MPI_ierr)
  int_tmp = MPI_orig(2)
  CALL MPI_GATHER(int_tmp,1,MPI_INT,MPI_block(1,3),1,MPI_INT,0,MPI_COMM_WORLD,MPI_ierr)
  int_tmp = MPI_orig(2)+ny
  CALL MPI_GATHER(int_tmp,1,MPI_INT,MPI_block(1,4),1,MPI_INT,0,MPI_COMM_WORLD,MPI_ierr)
  int_tmp = MPI_orig(3)
  CALL MPI_GATHER(int_tmp,1,MPI_INT,MPI_block(1,5),1,MPI_INT,0,MPI_COMM_WORLD,MPI_ierr)
  int_tmp = MPI_orig(3)+nz
  CALL MPI_GATHER(int_tmp,1,MPI_INT,MPI_block(1,6),1,MPI_INT,0,MPI_COMM_WORLD,MPI_ierr)

END SUBROUTINE initDomain


SUBROUTINE finalDomain
IMPLICIT NONE

  IF (ALLOCATED(  pm))      DEALLOCATE(  pm)
  IF (ALLOCATED(MPI_block)) DEALLOCATE(MPI_block)
  WRITE(*,*) 'Deallocation success: domain_mod.'

END SUBROUTINE finalDomain


SUBROUTINE setSize
IMPLICIT NONE

INTEGER       :: mod_x,mod_y,mod_z

  ! Determine nx and ny
  nx = nx_t / MPI_dims(1)
  ny = ny_t / MPI_dims(2)
  nz = nz_t / MPI_dims(3)
  mod_x = MOD(nx_t,MPI_dims(1))
  mod_y = MOD(ny_t,MPI_dims(2))
  mod_z = MOD(nz_t,MPI_dims(3))
  MPI_orig(1) = MPI_cdnt(1)*nx
  MPI_orig(2) = MPI_cdnt(2)*ny
  MPI_orig(3) = MPI_cdnt(3)*nz
  IF (mod_x.NE.0) THEN
    IF (MPI_cdnt(1).GE.(MPI_dims(1)-mod_x)) THEN
      nx = nx + 1
      MPI_orig(1) = MPI_orig(1) + (MPI_cdnt(1)-MPI_dims(1)+mod_x)
    END IF
  END IF
  IF (mod_y.NE.0) THEN
    IF (MPI_cdnt(2).GE.(MPI_dims(2)-mod_y)) THEN
      ny = ny + 1
      MPI_orig(2) = MPI_orig(2) + (MPI_cdnt(2)-MPI_dims(2)+mod_y)
    END IF
  END IF
  IF (mod_z.NE.0) THEN
    IF (MPI_cdnt(3).GE.(MPI_dims(3)-mod_z)) THEN
      nz = nz + 1
      MPI_orig(3) = MPI_orig(3) + (MPI_cdnt(3)-MPI_dims(3)+mod_z)
    END IF
  END IF

END SUBROUTINE setSize


SUBROUTINE readPm
IMPLICIT NONE

INTEGER(ST), DIMENSION(:), ALLOCATABLE :: dummy

  ! Read in bulk solid
  IF (pm_file.NE.'NO') THEN
    OPEN (UNIT=41,FILE=pm_file,STATUS='unknown')

    IF (MPI_orig(3).GT.0) THEN
    DO k=1,MPI_orig(3)
    DO j=1,ny_t
      READ(41,*)
    END DO
    END DO
    END IF

  DO k=1,nz
    IF (MPI_orig(2).GT.0) THEN
    DO j=1,MPI_orig(2)
      READ(41,*)
    END DO
    END IF
  
    IF (MPI_orig(1).GT.0) THEN
      IF (.NOT.ALLOCATED(dummy)) ALLOCATE(dummy(1:MPI_orig(1)),STAT=MPI_ierr)
      DO j=1,ny
        READ(41,*) dummy,pm(1:nx,j,k)
      END DO
      IF (ALLOCATED(dummy)) DEALLOCATE(dummy)
    ELSE
      DO j=1,ny
        READ(41,*) pm(1:nx,j,k)
      END DO
    END IF

    IF ((MPI_orig(2)+ny).LT.ny_t) THEN
    DO j=MPI_orig(2)+ny+1,ny_t
      READ(41,*)
    END DO
    END IF
  END DO
  END IF

  ! Exchange bound value
  CALL pEXCG(pm)

  ! Correct boundary
  IF (MPI_cdnt(1).NE.0) bound(1)=0
  IF (MPI_cdnt(2).NE.0) bound(2)=0
  IF (MPI_cdnt(3).NE.0) bound(3)=0
  IF (MPI_cdnt(1).NE.(MPI_dims(1)-1)) bound(4)=0
  IF (MPI_cdnt(2).NE.(MPI_dims(2)-1)) bound(5)=0
  IF (MPI_cdnt(3).NE.(MPI_dims(3)-1)) bound(6)=0

  IF (bound(1).NE.0) pm(0,:,:) = pm(1,:,:)
  IF (bound(2).NE.0) pm(:,0,:) = pm(:,1,:)
  IF (bound(3).NE.0) pm(:,:,0) = pm(:,:,1)
  IF (bound(4).NE.0) pm(nx+1,:,:) = pm(nx,:,:)
  IF (bound(5).NE.0) pm(:,ny+1,:) = pm(:,ny,:)
  IF (bound(6).NE.0) pm(:,:,nz+1) = pm(:,:,nz)

  IF (bound(1).EQ.1) pm(0,:,:) = pm_solid
  IF (bound(2).EQ.1) pm(:,0,:) = pm_solid
  IF (bound(3).EQ.1) pm(:,:,0) = pm_solid
  IF (bound(4).EQ.1) pm(nx+1,:,:) = pm_solid
  IF (bound(5).EQ.1) pm(:,ny+1,:) = pm_solid
  IF (bound(6).EQ.1) pm(:,:,nz+1) = pm_solid

END SUBROUTINE readPm


SUBROUTINE printDomain
IMPLICIT NONE

  IF (MPI_id.EQ.0) THEN
  10 FORMAT (A)
  11 FORMAT (A,I2)
  12 FORMAT (A,3I8)
  13 FORMAT (A,I8)
  14 FORMAT (A, ES13.5E3)
  15 FORMAT (A,3ES13.5E3)
  16 FORMAT (I3,A,6I7)
  WRITE(*,10) '********** Initialize domain **********'
  WRITE(*,11) ' Domain dimesion:',nd
  WRITE(*,12) ' Domain size:    ',nx_t,ny_t,nz_t
  WRITE(*,15) ' Domain length:  ',lx_t,ly_t,lz_t
  WRITE(*,13) '   Time size:    ',nt
  WRITE(*,14) '   Time length:  ',lt
  WRITE(*,14) '   Time step:    ',dt 
  WRITE(*,14) 'Lattice length:  ',dx
  WRITE(*,14) 'Lattice speed:   ',uc
  WRITE(*,10) '*** Successfully initialize domain! ***'
  WRITE(*,*)

  WRITE(*,10) ' ID |  x_beg  x_end  y_beg  y_end  z_beg  z_end'
  WRITE(*,10) '-----------------------------------------------'
  DO iip=1,MPI_np
  WRITE(*,16) iip-1,' |',MPI_block(iip,:)
  END DO
  WRITE(*,10) '-----------------------------------------------'
  WRITE(*,*)

  END IF

END SUBROUTINE printDomain

END MODULE domain_mod
