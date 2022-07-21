!!!--- Module file ---!!!
!!!-Description:
!!!   LBM model information.
!!!-Details
!!!   Parameters
!!!     Discrete velocity
!!!     MRT matrix
!!!     Weight parameters
!!!   Subroutines
!!!     Compute macros (locally)
!!!     Compute equilibrium (locally)
!!!     Perform collision (locally)
!!!-Date
!!!   2018-11-27   Eric      V2.0
!!!-Declaration
!!!   LBM model module


MODULE lbmodel_mod
USE    mpi_mod, only: MPI_id
USE  const_mod, only: nd,ST,DBL
IMPLICIT NONE

!*** Velocity dimension ***!
INTEGER(ST), PARAMETER       :: nq=15  ! degree of freedom (DOF) of discrete velocity
INTEGER(ST)                  :: iiq    ! lattice velocity index


!*** Discrete vilocity ***!
INTEGER(ST), DIMENSION(nq,nd):: e      ! discrete lattice velocity
INTEGER(ST), DIMENSION(nq)   :: opp    ! oppsoite velocity index
REAL(DBL)                    :: RT     ! p=rho * RT, RT=cs**2, here RT=(cs/c)**2


!*** MRT ***!
REAL(DBL), DIMENSION(nq,nq)  :: M,M_inv     ! MRT matrix and its inv
REAL(DBL), DIMENSION(nq,nq)  :: K_diag      ! diagonal matrix of relxation factor
REAL(DBL)                    :: w_tau       ! reciprocal of relaxation time


!*** Weight parameters ***!
! Equilibrium distribution: wt
!   f=rho*wt*(1+(u and ci))
REAL(DBL), DIMENSION(nq)     :: wt


CONTAINS


SUBROUTINE initLBModel
IMPLICIT NONE

CHARACTER(70) :: dummy

  ! Input
  OPEN (UNIT=41,FILE='INPUT/lbmodelDict',STATUS='UNKNOWN')
  READ(41,*)
  READ(41,*) dummy,w_tau
  CLOSE (UNIT=41)

  e(:,1) = (/1,0,0,-1, 0, 0,1,-1, 1,-1, 1,-1, 1,-1,0/)
  e(:,2) = (/0,1,0, 0,-1, 0,1, 1,-1,-1, 1, 1,-1,-1,0/)
  e(:,3) = (/0,0,1, 0, 0,-1,1, 1, 1, 1,-1,-1,-1,-1,0/)

  opp = (/4,5,6,1,2,3,14,13,12,11,10,9,8,7,15/)

  RT = 1d0/3d0

  ! D3Q15
  ! M(1:3 ,:) are related with momentum rho*ux,rho*uy
  ! M(6:10,:) are related with viscosity
  M(15,:) = (/ 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0/) ! 1  - rho - 0  0
  M( 1,:) = (/ 1d0, 0d0, 0d0,-1d0, 0d0, 0d0, 1d0,-1d0, 1d0,-1d0, 1d0,-1d0, 1d0,-1d0, 0d0/) ! ex - ux  - 0  0
  M( 2,:) = (/ 0d0, 1d0, 0d0, 0d0,-1d0, 0d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 0d0/) ! ey - uy  - 0  0
  M( 3,:) = (/ 0d0, 0d0, 1d0, 0d0, 0d0,-1d0, 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0, 0d0/) ! ez - uz  - 0  0
  M( 4,:) = (/-1d0,-1d0,-1d0,-1d0,-1d0,-1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0,-2d0/) !          - se s1
  M( 5,:) = (/-4d0,-4d0,-4d0,-4d0,-4d0,-4d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0,16d0/) !          - sn s2
  M( 6,:) = (/ 2d0,-1d0,-1d0, 2d0,-1d0,-1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/) !          - sv s9
  M( 7,:) = (/ 0d0, 1d0,-1d0, 0d0, 1d0,-1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/) !          - sv s9
  M( 8,:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 1d0, 0d0/) !          - sv s11
  M( 9,:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0, 1d0, 1d0, 0d0/) !          - sv s11
  M(10,:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0,-1d0, 1d0,-1d0,-1d0, 1d0,-1d0, 1d0, 0d0/) !          - sv s11
  M(11,:) = (/-4d0, 0d0, 0d0, 4d0, 0d0, 0d0, 1d0,-1d0, 1d0,-1d0, 1d0,-1d0, 1d0,-1d0, 0d0/) !          - sq s4
  M(12,:) = (/ 0d0,-4d0, 0d0, 0d0, 4d0, 0d0, 1d0, 1d0,-1d0,-1d0, 1d0, 1d0,-1d0,-1d0, 0d0/) !          - sq s4
  M(13,:) = (/ 0d0, 0d0,-4d0, 0d0, 0d0, 4d0, 1d0, 1d0, 1d0, 1d0,-1d0,-1d0,-1d0,-1d0, 0d0/) !          - sq s4
  M(14,:) = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0,-1d0,-1d0, 1d0,-1d0, 1d0, 1d0,-1d0, 0d0/) !          - sm s14

  DO iiq=1,nq
    M_inv(iiq,:) = M(iiq,:) / DOT_PRODUCT(M(iiq,:),M(iiq,:))
  END DO
  M_inv = transpose(M_inv)

  CALL initMRT

  ! D3Q15
  wt = (/ 8d0, 8d0, 8d0, 8d0, 8d0, 8d0, 1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0, 16d0/) / 72d0

END SUBROUTINE initLBModel


SUBROUTINE initMRT
IMPLICIT NONE

  ! Diagonal matrix - related to bulk properties
  !--This MRT setting is from paper:
  !--  Chongxun Pan-Computers & Fluids 35 (2006)
  !--Set A:
  !--  sv = w_tau; se=sn=sp=sm=sq=8*(2-sv)/(8-sv)
  !--Set B:
  !--  sv=se=sn=sp=w_tau; sm=sq=8*(2-sv)/(8-sv)
  !--Present is set A
  K_diag = 0d0
  DO iiq=4,nq-1
    K_diag(iiq,iiq) = 8d0*(2d0-w_tau)/(8d0-w_tau)
  END DO
  DO iiq =6,10
    K_diag(iiq,iiq) = w_tau
  END DO
!  !--Presnet is Set B
!  K_diag = 0d0
!  DO iiq=4,nq-1
!    K_diag(iiq,iiq) = w_tau
!  END DO
!  DO iiq =11,14
!    K_diag(iiq,iiq) = 8d0*(2d0-w_tau)/(8d0-w_tau)
!  END DO

END SUBROUTINE initMRT

SUBROUTINE printLBModel
IMPLICIT NONE
INTEGER(ST)   :: jjq,iid

  ! Print lbm model information
  IF (MPI_id.EQ.0) THEN
  10 FORMAT (A)
  11 FORMAT (A,I2)
  12 FORMAT (15I3)
  13 FORMAT (15F7.3)
  14 FORMAT (A,15F7.3)
  15 FORMAT (A,F7.3)
  WRITE(*,10) '********* Initailize LBM model ********'
  WRITE(*,11) 'Number of discrete velocity:',nq
  WRITE(*,10) 'Discrete velocity:'
  WRITE(*,12) ((e(iiq,iid),iiq=1,nq),iid=1,nd)
  WRITE(*,10) 'MRT matrix:'
  WRITE(*,13) ((M(iiq,jjq),jjq=1,nq),iiq=1,nq)
  WRITE(*,10) 'MRT inv matrix:'
  WRITE(*,13) ((M_inv(iiq,jjq),jjq=1,nq),iiq=1,nq)
  WRITE(*,10) 'Weight parameters:'
  WRITE(*,14) 'wt: ',(wt(iiq),iiq=1,nq)
  WRITE(*,15) 'Relaxation time:               ',1d0/w_tau
  !WRITE(*,15) 'Scaling parameter of K diag  : ',wK
  WRITE(*,10) '** Successfully initailize LBM model **'
  WRITE(*,*)
  END IF

END SUBROUTINE printLBModel


SUBROUTINE collision(f1,force_bd)
IMPLICIT NONE
    
REAL(DBL), INTENT(INOUT),DIMENSION(nq)      :: f1
REAL(DBL), INTENT(IN   ),DIMENSION(nd)      :: force_bd

REAL(DBL), DIMENSION(nq)      :: feq1
REAL(DBL), DIMENSION(nd)      :: u1
REAL(DBL)                     :: rho1

  CALL computeRhoU(rho1,u1,f1)
  CALL computeEqlm(feq1,rho1,u1)
  CALL collisionSGP(f1,feq1,force_bd)

END SUBROUTINE collision


!*** collisionSGP ***!
!-Description
!   Single phase collision (locally)
!-Arguments
!   f1:   distribution
!   feq1: equilibrium distribution
!   force_bd:  body force
SUBROUTINE collisionSGP(f1,feq1,force_bd)
USE  domain_mod, only: iit
IMPLICIT NONE
    
REAL(DBL), INTENT(INOUT),DIMENSION(nq)      :: f1
REAL(DBL), INTENT(IN),   DIMENSION(nq)      :: feq1
REAL(DBL), INTENT(IN),   DIMENSION(nd), OPTIONAL :: force_bd

REAL(DBL), DIMENSION(nq)     :: f_src  ! source term related with body force

  f_src = 0d0
  IF (iit.GE.0) f_src(1:3) = SUM(f1)*force_bd(:)

  f1 = f1 - MATMUL(M_inv,MATMUL(K_diag,MATMUL(M,f1-feq1))) + MATMUL(M_inv,f_src)

END SUBROUTINE collisionSGP


!*** computeEqlm ***!
!-Description
!   Compute equilibrium distribution (locally)
!-Arguments
!   feq1: equilibrium distribution
!   rho1: density
!   u1:   velocity
SUBROUTINE computeEqlm(feq1,rho1,u1)
IMPLICIT NONE

REAL(DBL), INTENT(OUT), DIMENSION(nq)  :: feq1
REAL(DBL), INTENT(IN)                  :: rho1
REAL(DBL), INTENT(IN),  DIMENSION(nd)  :: u1

REAL(DBL)     :: e_dot_u, u_sqr, ue_temp

  u_sqr = DOT_PRODUCT(u1,u1)      ! square velocity

  feq1 = 0d0
  DO iiq=1,nq
    e_dot_u    = DOT_PRODUCT(e(iiq,:),u1)        ! dot production of e and u
    feq1(iiq)  = rho1*wt(iiq)*(1d0 + 3d0*e_dot_u + 4.5d0*e_dot_u**2 - 1.5d0*u_sqr)
  END DO

END SUBROUTINE computeEqlm


!*** computeRho ***!
!-Description
!   Compute density (locally).
!-Arguments
!  rho1: density
!  f1:   distribution
SUBROUTINE computeRho(rho1,f1)
IMPLICIT NONE

REAL(DBL), INTENT(OUT)                 :: rho1
REAL(DBL), INTENT(IN), DIMENSION(nq)   :: f1

  rho1 = SUM(f1)

END SUBROUTINE computeRho


!*** computeU ***!
!-Description
!   Compute velocity (locally).
!-Arguments
!   u1:   velocity
!   rho1: density
!   f1:   distribution
SUBROUTINE computeU(u1,rho1,f1)
IMPLICIT NONE

REAL(DBL), INTENT(OUT),DIMENSION(nd)   :: u1
REAL(DBL), INTENT(IN)                  :: rho1
REAL(DBL), INTENT(IN), DIMENSION(nq)   :: f1

  u1 = MATMUL(f1,e) / rho1

END SUBROUTINE computeU


!*** computeRhoU ***!
!-Description
!   Compute density and velocity (locally).
!-Arguments
!   rho1: density
!   u1:   velocity
!   f1:   distribution
SUBROUTINE computeRhoU(rho1,u1,f1)
IMPLICIT NONE

REAL(DBL), INTENT(OUT)                 :: rho1
REAL(DBL), INTENT(OUT),DIMENSION(nd)   :: u1
REAL(DBL), INTENT(IN), DIMENSION(nq)   :: f1

  rho1 = SUM(f1)
  u1 = MATMUL(f1,e) / rho1

END SUBROUTINE computeRhoU

END MODULE lbmodel_mod
