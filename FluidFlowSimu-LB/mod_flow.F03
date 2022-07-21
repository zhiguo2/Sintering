!!!--- Module file ---!!!
!!!-Description:
!!!   Flow fields information.
!!!-Details
!!!   Fluid properties
!!!   Field quantities
!!!   Body forces
!!!   Pressure / velocity boundary value
!!!-Date
!!!   2018-11-27   Eric      V2.0
!!!-Declaration
!!!   Fields / Main operation data module


MODULE flow_mod
USE     mpi_mod, only: MPI_id
USE   const_mod, only: nd,ST,DBL,PI
USE  domain_mod, only: nx,ny,nz
USE lbmodel_mod, only: nq
IMPLICIT NONE

!*** Fluid properties ***!
REAL(DBL)     :: rho0   ! [kg/m3] density
REAL(DBL)     :: nu0    ! [m2/s] viscosity


!*** Fields ***!
REAL(DBL), DIMENSION(:,:,:,:), ALLOCATABLE  :: f      ! distribution function
REAL(DBL), DIMENSION(:,:,:,:), ALLOCATABLE  :: u      ! velocity
REAL(DBL), DIMENSION(  :,:,:), ALLOCATABLE  :: rho    ! density


!*** Body force ***!
REAL(DBL), DIMENSION(nd)     :: force_bd    ! body force


!*** Boundary ***!
INTEGER, PARAMETER :: bI_cyclic         = 0
INTEGER, PARAMETER :: bI_bounce_back    = 1
INTEGER, PARAMETER :: bI_symmetry       = 2
INTEGER, PARAMETER :: bI_fixed_pressure = 3
INTEGER, PARAMETER :: bI_inlet_velocity = 4
REAL(DBL), DIMENSION(6)      ::  um_bound = 0d0  ! inlet velocity
REAL(DBL), DIMENSION(6)      ::   p_bound = 0d0  ! in/outlet pressure


CONTAINS


SUBROUTINE initFlow
USE   const_mod, only: nd,DBL,PI,RGC
USE     mpi_mod, only: MPI_ierr,subdir
USE  domain_mod, only: i,j,k,dx,dy,dz,dt,nt,lt,uc,nt_s,MPI_orig
USE lbmodel_mod, only: RT,w_tau,computeEqlm
IMPLICIT NONE

CHARACTER(90) :: dummy
CHARACTER(7)  :: fileI
INTEGER       :: init_from_file
REAL(DBL)     :: p_ir,p_or
REAL(DBL),DIMENSION(nd) :: para0,para1

  ! Input
  OPEN (UNIT=41,FILE='INPUT/flowDict',STATUS='UNKNOWN')
  READ(41,*)
  READ(41,*)
  READ(41,*) dummy,rho0
  READ(41,*) dummy,nu0  ! read in kinetic viscosity
  READ(41,*) dummy,force_bd
  READ(41,*)
  READ(41,*)
  READ(41,*) dummy,  p_bound
  READ(41,*) dummy, um_bound
  READ(41,*)
  READ(41,*)
  READ(41,*) dummy,init_from_file 
  IF (init_from_file.NE.0) THEN
  READ(41,*) dummy,para0
  READ(41,*) dummy,para1
  READ(41,*) dummy,p_ir,p_or
  END IF
  CLOSE (UNIT=41)

  dt = dx**2/nu0/3d0*(1d0/w_tau-0.5d0)
  uc = dx/dt
  lt = nt*dt

  um_bound = um_bound /  uc
   p_bound =  p_bound /  uc**2
  !     nu0 =      nu0 / (uc*dx)
  force_bd = force_bd / (uc**2/dx)

  ! Allocate memory
  IF (.NOT.ALLOCATED(    f)) ALLOCATE(    f(nq,0:(nx+1),0:(ny+1),0:(nz+1)),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(    u)) ALLOCATE(    u(nd,0:(nx+1),0:(ny+1),0:(nz+1)),STAT=MPI_ierr)
  IF (.NOT.ALLOCATED(  rho)) ALLOCATE(  rho(   0:(nx+1),0:(ny+1),0:(nz+1)),STAT=MPI_ierr)

  ! Initialize field
  IF (init_from_file.EQ.0) THEN    

    WRITE(fileI,'(I7.7)') nt_s
    dummy = 'OUTPUT/'//TRIM(subdir(MPI_id))//'/'//TRIM(fileI)//'/f.dat'
    OPEN (UNIT=40, FILE=TRIM(dummy),STATUS='UNKNOWN',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
    DO k=0,nz+1
    DO j=0,ny+1
    DO i=0,nx+1
      READ(40) f(:,i,j,k)
    END DO
    END DO
    END DO
    CLOSE(UNIT=40)

  ELSE

    ! Initialize color function
    rho = rho0 + p_or/RT/uc**2
    IF (init_from_file.EQ.1) THEN
      DO k=0,nz+1
      DO j=0,ny+1
      DO i=0,nx+1
        IF ((((MPI_orig(1)+i-0.5d0-para0(1)/dx)/para1(1)*dx)**2+&
            &((MPI_orig(2)+j-0.5d0-para0(2)/dy)/para1(2)*dy)**2+&
            &((MPI_orig(3)+k-0.5d0-para0(3)/dz)/para1(3)*dz)**2).LE.1d0) THEN
          rho(i,j,k) = rho0 + p_ir/RT/uc**2
        END IF
      END DO
      END DO
      END DO
    ELSE IF (init_from_file.EQ.2) THEN
      DO k=0,nz+1
      DO j=0,ny+1
      DO i=0,nx+1
        IF (((MPI_orig(1)+i-0.5d0).GE.(para0(1)/dx)).AND.&
           &((MPI_orig(2)+j-0.5d0).GE.(para0(2)/dy)).AND.&
           &((MPI_orig(3)+k-0.5d0).GE.(para0(3)/dz)).AND.&
           &((MPI_orig(1)+i-0.5d0).LE.(para1(1)/dx)).AND.&
           &((MPI_orig(2)+j-0.5d0).LE.(para1(2)/dy)).AND.&
           &((MPI_orig(3)+k-0.5d0).LE.(para1(3)/dz))) THEN
          rho(i,j,k) = rho0 + p_ir/RT/uc**2
        END IF
      END DO
      END DO
      END DO
    END IF
  
    ! Initialize u
    u = 0d0

    ! Initialize distribution function with equilibrium distribution
    DO k=0,nz+1
    DO j=0,ny+1
    DO i=0,nx+1
      CALL computeEqlm(f(:,i,j,k),rho(i,j,k),u(:,i,j,k))
    END DO
    END DO
    END DO

  END IF

END SUBROUTINE initFlow


SUBROUTINE printFlow
IMPLICIT NONE

  IF (MPI_id.EQ.0) THEN
  10 FORMAT (A)
  11 FORMAT (A,I2)
  12 FORMAT (A, ES14.5E3)
  13 FORMAT (A,2ES14.5E3)
  14 FORMAT (A,3ES14.5E3)
  WRITE(*,10) '*********** Initialize fluid **********'
  WRITE(*,13) '  Density: ',rho0
  WRITE(*,13) 'Viscosity: ', nu0
  WRITE(*,14) 'Body force:',force_bd
  WRITE(*,10) '**** Succeddfully initialize fluid ****'
  WRITE(*,*)
  END IF

END SUBROUTINE printFlow


SUBROUTINE finalFlow
IMPLICIT NONE

  IF (ALLOCATED(    f)) DEALLOCATE(    f)
  IF (ALLOCATED(    u)) DEALLOCATE(    u)
  IF (ALLOCATED(  rho)) DEALLOCATE(  rho)

  WRITE(*,*) 'Deallocate success: flow_mod!'

END SUBROUTINE finalFlow


SUBROUTINE writeFieldInfo
USE    mpi_mod, only: pMAXVAL,pMINVAL
USE domain_mod, only: dx,dt,uc

REAL(DBL)     :: min_f,min_rho,min_u
REAL(DBL)     :: max_f,max_rho,max_u

  ASSOCIATE (    f_f =>     f(:,1:nx,1:ny,1:nz),&
            &  rho_f =>   rho(  1:nx,1:ny,1:nz),&
            &    u_f =>     u(:,1:nx,1:ny,1:nz))
  min_f   = pMINVAL(f_f)
  max_f   = pMAXVAL(f_f)
  min_rho = pMINVAL(rho_f)
  max_rho = pMAXVAL(rho_f)
  min_u   = pMINVAL(u_f)*uc
  max_u   = pMAXVAL(u_f)*uc

  ! Print field information
  IF (MPI_id.EQ.0) THEN
  10 FORMAT (A)
  12 FORMAT (A,ES13.5E3,1X,A,ES13.5E3)
  13 FORMAT (A,ES13.5E3)
  WRITE(*,10) '-----------------------------------------'
  WRITE(*,10) ' Field name |   Min value  |   Max value'
  WRITE(*,10) '-----------------------------------------'
  WRITE(*,12) '       f    |', min_f    ,'|', max_f
  WRITE(*,12) '     rho    |', min_rho  ,'|', max_rho
  WRITE(*,12) '       u    |', min_u    ,'|', max_u
  WRITE(*,10) '-----------------------------------------'
  WRITE(*,*)
  END IF
  END ASSOCIATE

END SUBROUTINE writeFieldInfo


!*** stream ***!
!-Description
!   Perform stream process globally
!-Arguments
!   f: distribution function
SUBROUTINE stream(f)
USE lbmodel_mod, only: e,iiq
IMPLICIT NONE

REAL(DBL), INTENT(INOUT), DIMENSION(nq,0:(nx+1),0:(ny+1),0:(nz+1))  :: f

  DO iiq=1,nq-1
    f(iiq,:,:,:) = CSHIFT(f(iiq,:,:,:),-e(iiq,1),1)
    f(iiq,:,:,:) = CSHIFT(f(iiq,:,:,:),-e(iiq,2),2)
    f(iiq,:,:,:) = CSHIFT(f(iiq,:,:,:),-e(iiq,3),3)
  END DO
  
END SUBROUTINE stream

END MODULE flow_mod
