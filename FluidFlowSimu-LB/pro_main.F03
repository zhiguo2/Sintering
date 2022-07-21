!!!--- Program file ---!!!
!!!-Description:
!!!   MPI program for 3D onephase flow process.
!!!-Details
!!!   FLOWLBM3D
!!!     Flow LBM program in 3D.
!!!-Date
!!!   2018-11-27   Eric      V2.0
!!!-Declaration
!!!   Main program


!*** Main program ***!
!-Description
PROGRAM FLOWLBM3D
USE   const_mod
USE     mpi_mod
USE  domain_mod
USE lbmodel_mod
USE    flow_mod
USE   bound_mod
USE      io_mod
IMPLICIT NONE

REAL(DBL)     :: time_beg,time_end
REAL(DBL), DIMENSION(nq)     :: feq1

CALL CPU_TIME(time_beg)

CALL initCase

IF (MPI_id.EQ.0) THEN
  WRITE(*,*) '**************** Start computation ****************'
  WRITE(*,*)
END IF

!*** main loop ***!
DO iit = nt_s,nt_s+nt

  IF (MPI_id.EQ.0) THEN
    WRITE(*,'(A,I9)') 'Time step: ',iit
  END IF

  ! Compute macro properties for each phase
  DO k=1,nz
  DO j=1,ny
  DO i=1,nx
  IF (pm(i,j,k).EQ.pm_pore) THEN
    CALL computeRhoU(rho(i,j,k),u(:,i,j,k),   f(:,i,j,k))
  END IF
  END DO
  END DO
  END DO

  ! Collison
  DO k=1,nz
  DO j=1,ny
  DO i=1,nx
  IF (pm(i,j,k).EQ.pm_pore) THEN
    ! compute equilibrium distribution
    CALL computeEqlm(feq1 ,rho(i,j,k),    u(:,i,j,k))
    !IF (MINVAL(feq1).LT.0) WRITE(*,*) 'Negative equilibrium distribution !!!'
    
    ! Single-phase collision
    CALL collisionSGP(f(:,i,j,k),feq1,force_bd)
    !IF (MINVAL(f(:,i,j,k)).LT.0) WRITE(*,*) 'Negative distribution after collision !!!'
  END IF 
  END DO
  END DO
  END DO

  ! update ghost layer
  !-- This update is mainly due to the boundary condition.
  !-- For fixed pressure boundary, f and rho in ghost layer may be needed.
  CALL pEXCG(rho)
  CALL pEXCG(u)
  CALL pEXCG(f)
  
  ! inlet/outlet boundary
  DO iib=1,6
    IF (bound(iib).EQ.bI_fixed_pressure) THEN
    CALL fixedPressure(p_bound(iib) ,iib)
    END IF

    IF (bound(iib).EQ.bI_inlet_velocity) THEN
    CALL inletVelocity(um_bound(iib) ,iib)
    END IF
  END DO

  ! symmetry boundary
  DO iib=1,6
    IF (bound(iib).EQ.bI_symmetry) THEN
    CALL symmetry(f,iib)
    END IF
  END DO

  ! Perform streaming of each phase
  CALL stream(f)
  CALL bounceBack(f)
  CALL bounceBackB(f)

  ! Output time series results
  IF ((MOD(iit-nt_s,nt_o).EQ.0) .AND. (iit.GE.0)) THEN
    CALL writeFieldInfo
    CALL outputVTK
  END IF

  IF ((MOD(iit-nt_s,nt_f).EQ.0) .AND. (iit.GT.0)) THEN
    CALL outputF
  END IF

  CALL outputFile

END DO

CALL finalCase

CALL MPI_FINALIZE(MPI_ierr)

CALL CPU_TIME(time_end)
WRITE(*,*) 'Elapsed time is ', time_end-time_beg, ' s.'

END PROGRAM FLOWLBM3D
