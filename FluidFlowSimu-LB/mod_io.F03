!!!--- Module file ---!!!
!!!-Description:
!!!   Subroutines for input and output.
!!!-Details
!!!   Subroutine
!!!     Initialize each module.
!!!     Finalize dynamic arrays.
!!!     Output vtk files.
!!!     Output intermediate results.
!!!-Date
!!!   2018-11-27   Eric      V2.0
!!!-Declaration
!!!   IO module


MODULE io_mod
IMPLICIT NONE


CONTAINS


!*** initializeCase ***!
!-Description
!   Initialize each module and fields.
!   Fields are intialized with external subroutine - initializeField.
!-Arguments
!   NO
SUBROUTINE initCase
USE     mpi_mod, only: initMPI,    printMPI,MPI_id
USE  domain_mod, only: initDomain, printDomain
USE lbmodel_mod, only: initLBModel,printLBModel
USE    flow_mod, only: initFlow,   printFlow
IMPLICIT NONE

  CALL initMPI

  IF (MPI_id.EQ.0) THEN
  WRITE(*,*)
  WRITE(*,*) '*********** Start initialization ***********'
  END IF

  CALL initDomain
  CALL initLBModel
  CALL initFlow

  CALL printDomain
  CALL printLBModel
  CALL printFlow

  IF (MPI_id.EQ.0) THEN
  WRITE(*,*) '*********** Finish initialization **********'
  WRITE(*,*)
  END IF

END SUBROUTINE initCase


SUBROUTINE finalCase
USE domain_mod, only: finalDomain
USE   flow_mod, only: finalFlow
IMPLICIT NONE

  CALL finalDomain
  CALL finalFlow

END SUBROUTINE finalCase


!*** outputVTK ***!
!-Description
!   Output time instant results into VTK files.
!-Arguments
!   NO
SUBROUTINE outputVTK
USE  domain_mod
USE     mpi_mod
USE lbmodel_mod, only: RT
USE    flow_mod
IMPLICIT NONE

INTEGER(LG)   :: x0,x1,y0,y1,z0,z1,n_byte
CHARACTER(90) :: dummy
CHARACTER(7)  :: fileI
CHARACTER     :: LF

  LF = ACHAR(10)

  x0 = MPI_orig(1)
  y0 = MPI_orig(2)
  z0 = MPI_orig(3)
  x1 = MPI_orig(1)+nx
  y1 = MPI_orig(2)+ny
  z1 = MPI_orig(3)+nz

  WRITE(fileI,'(I7.7)') iit

  ! Output to serial VTK XML file with each processor's data
  dummy = 'OUTPUT/'//TRIM(subdir(MPI_id))//'/result.'//TRIM(fileI)//'.vtr'
  OPEN (UNIT=40,FILE=TRIM(dummy),STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
  WRITE(40) '<?xml version="1.0"?>'//LF
  WRITE(40) '<VTKFile type="RectilinearGrid" Version="1.0"'//LF
  WRITE(40) '         byte_order="BigEndian" header_type="UInt64">'//LF
  WRITE(dummy,'(I0,X,I0,X,I0,X,I0,X,I0,X,I0)') x0,x1,y0,y1,z0,z1
  WRITE(40) '  <RectilinearGrid WholeExtent="'//TRIM(dummy)//'">'//LF
  WRITE(40) '    <Piece Extent="'//TRIM(dummy)//'">'//LF
  WRITE(40) '      <PointData>'//LF
  WRITE(40) '      </PointData>'//LF
  WRITE(40) '      <CellData Scalars="p s" Vectors="u">'//LF
  WRITE(40) '        <DataArray type="Float64" Name="p" format="appended" offset="0">'//LF
  WRITE(40) '        </DataArray>'//LF
  n_byte = SIZEOF(REAL(1,DBL))*nx*ny*nz + 8
  WRITE(dummy,'(I0)') n_byte
  WRITE(40) '        <DataArray type="Int8" Name="s" format="appended" offset="'//TRIM(dummy)//'">'//LF
  WRITE(40) '        </DataArray>'//LF
  n_byte = SIZEOF(INT(1,ST))*nx*ny*nz + 8 + n_byte
  WRITE(dummy,'(I0)') n_byte
  WRITE(40) '        <DataArray type="Float64" Name="u" NumberOfComponents="3" format="appended"'//LF
  WRITE(40) 'offset="'//TRIM(dummy)//'">'//LF
  WRITE(40) '        </DataArray>'//LF
  WRITE(40) '      </CellData>'//LF
  WRITE(40) '      <Coordinates>'//LF
  n_byte = SIZEOF(REAL(1,DBL))*nx*ny*nz*3 + 8 + n_byte
  WRITE(dummy,'(I0)') n_byte
  WRITE(40) '        <DataArray type="Float64" Name="x" format="appended" offset="'//TRIM(dummy)//'">'//LF
  WRITE(40) '        </DataArray>'//LF
  n_byte = SIZEOF(REAL(1,DBL))*(nx+1) + 8 + n_byte
  WRITE(dummy,'(I0)') n_byte
  WRITE(40) '        <DataArray type="Float64" Name="y" format="appended" offset="'//TRIM(dummy)//'">'//LF
  WRITE(40) '        </DataArray>'//LF
  n_byte = SIZEOF(REAL(1,DBL))*(ny+1) + 8 + n_byte
  WRITE(dummy,'(I0)') n_byte
  WRITE(40) '        <DataArray type="Float64" Name="z" format="appended" offset="'//TRIM(dummy)//'">'//LF
  WRITE(40) '        </DataArray>'//LF
  WRITE(40) '      </Coordinates>'//LF
  WRITE(40) '    </Piece>'//LF
  WRITE(40) '  </RectilinearGrid>'//LF
  WRITE(40) '  <AppendedData encoding="raw">'//LF
  WRITE(40) '    _'
  n_byte = SIZEOF(REAL(1,DBL))*nx*ny*nz
  WRITE(40) n_byte,(((RT*(rho(i,j,k)-rho0)*uc**2,i=1,nx),j=1,ny),k=1,nz)
  n_byte = SIZEOF(INT(1,ST))*nx*ny*nz
  WRITE(40) n_byte,(((pm(i,j,k),i=1,nx),j=1,ny),k=1,nz)
  n_byte = SIZEOF(REAL(1,DBL))*nx*ny*nz*3
  WRITE(40) n_byte,(((u(:,i,j,k)*uc,i=1,nx),j=1,ny),k=1,nz)
  n_byte = SIZEOF(REAL(1,DBL))*(nx+1)
  WRITE(40) n_byte,((MPI_orig(1)+i)*dx,i=0,nx)
  n_byte = SIZEOF(REAL(1,DBL))*(ny+1)
  WRITE(40) n_byte,((MPI_orig(2)+j)*dy,j=0,ny)
  n_byte = SIZEOF(REAL(1,DBL))*(nz+1)
  WRITE(40) n_byte,((MPI_orig(3)+k)*dz,k=0,nz)
  WRITE(40) LF//'  </AppendedData>'//LF
  WRITE(40) '</VTKFile>'//LF
  
  CLOSE(UNIT=40)   ! close file

  IF (MPI_id.EQ.0) THEN
  ! Output to parallel VTK XML file with partition info
  OPEN (UNIT=41, FILE='OUTPUT/result.'//trim(fileI)//'.pvtr',STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
  WRITE(41) '<?xml version="1.0"?>'//LF
  WRITE(41) '<VTKFile type="PRectilinearGrid" Version="1.0"'//LF
  WRITE(41) '         byte_order="BigEndian"  header_type="UInt64">'//LF
  WRITE(dummy,'(I0,X,I0,X,I0,X,I0,X,I0,X,I0)') 0,nx_t,0,ny_t,0,nz_t
  WRITE(41) '  <PRectilinearGrid WholeExtent="'//TRIM(dummy)//'"'//LF
  WRITE(41) '                    GhostLevel="0">'//LF
  WRITE(41) '    <PPointData>'//LF
  WRITE(41) '    </PPointData>'//LF
  WRITE(41) '    <PCellData Scalars="p s" Vectors="u">'//LF
  WRITE(41) '      <PDataArray type="Float64" Name="p"/>'//LF
  WRITE(41) '      <PDataArray type="Int8" Name="s"/>'//LF
  WRITE(41) '      <PDataArray type="Float64" Name="u" NumberOfComponents="3"/>'//LF
  WRITE(41) '    </PCellData>'//LF
  WRITE(41) '    <PCoordinates>'//LF
  WRITE(41) '      <PDataArray type="Float64" Name="x"/>'//LF
  WRITE(41) '      <PDataArray type="Float64" Name="y"/>'//LF
  WRITE(41) '      <PDataArray type="Float64" Name="z"/>'//LF
  WRITE(41) '    </PCoordinates>'//LF
  DO iip=1,MPI_np
  WRITE(dummy,'(I0,X,I0,X,I0,X,I0,X,I0,X,I0)') MPI_block(iip,:)
  WRITE(41) '    <Piece Extent="'//TRIM(dummy)//'"'//LF
  WRITE(41) '           Source="'//TRIM(subdir(iip-1))//'/result.'//TRIM(fileI)//'.vtr"/>'//LF
  END DO
  WRITE(41) '  </PRectilinearGrid>'//LF
  WRITE(41) '</VTKFile>'//LF
  
  CLOSE(UNIT=41)   ! close file
  END IF

END SUBROUTINE outputVTK


!*** outputF ***!
!-Description
!   Output distribution function results to dat file for continuous computation.
!-Arguments
!   NO
SUBROUTINE outputF
USE  domain_mod
USE     mpi_mod
USE    flow_mod
IMPLICIT NONE

CHARACTER(90) :: dummy
CHARACTER(7)  :: fileI
INTEGER       :: iiq

  WRITE(fileI,'(I7.7)') iit

  CALL SYSTEM('mkdir -p OUTPUT/'//TRIM(subdir(MPI_id))//'/'//TRIM(fileI))

  dummy = 'OUTPUT/'//TRIM(subdir(MPI_id))//'/'//TRIM(fileI)//'/f.dat'
  OPEN (UNIT=40, FILE=TRIM(dummy),STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
  WRITE(40) ((((f(iiq,i,j,k), iiq=1,nq), i=0,nx+1), j=0,ny+1), k=0,nz+1)
  CLOSE(UNIT=40)

END SUBROUTINE outputF


SUBROUTINE outputFile
USE   const_mod, only: DBL
USE     mpi_mod
USE lbmodel_mod, only: RT
USE  domain_mod, only: dx,dt,uc,nx,ny,nz,iit,pm,pm_pore
USE    flow_mod, only: f,rho,u,rho0
IMPLICIT NONE

REAL(DBL)     :: min_f,min_rho,min_p,min_u
REAL(DBL)     :: max_f,max_rho,max_p,max_u
REAL(DBL)     :: sum_ufx ,sum_ufy ,sum_ufz
REAL(DBL)     :: sum_ufxi,sum_ufyi,sum_ufzi
INTEGER       :: fUnit
LOGICAL       :: ok

  ASSOCIATE (    f_f =>   f(:,1:nx,1:ny,1:nz),&
            &  rho_f => rho(  1:nx,1:ny,1:nz),&
            &    u_f =>   u(:,1:nx,1:ny,1:nz),&
            &   pm_f =>  pm(  1:nx,1:ny,1:nz) )
  min_f   = pMINVAL(f_f)
  max_f   = pMAXVAL(f_f)
  min_rho   = pMINVAL(rho_f)
  max_rho = pMAXVAL(rho_f)
  min_u   = pMINVAL(u_f)*uc
  max_u   = pMAXVAL(u_f)*uc

  max_p = RT*(max_rho-rho0)*uc**2
  min_p = RT*(min_rho-rho0)*uc**2

  sum_ufxi = SUM(u_f(1,:,:,:),MASK=pm_f.EQ.pm_pore) *uc*dx**3
  sum_ufyi = SUM(u_f(2,:,:,:),MASK=pm_f.EQ.pm_pore) *uc*dx**3
  sum_ufzi = SUM(u_f(3,:,:,:),MASK=pm_f.EQ.pm_pore) *uc*dx**3
  CALL MPI_REDUCE(sum_ufxi,sum_ufx,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPI_ierr)
  CALL MPI_REDUCE(sum_ufyi,sum_ufy,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPI_ierr)
  CALL MPI_REDUCE(sum_ufzi,sum_ufz,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,MPI_ierr)
  END ASSOCIATE

  IF (MPI_id.EQ.0) THEN
  20 FORMAT (A)
  21 FORMAT (I7,9ES14.5E3)
  INQUIRE(FILE='fieldInfo.dat',OPENED=ok,NUMBER=fUnit)
  IF (ok) THEN
    WRITE(fUnit,21) iit,sum_ufx,sum_ufy,sum_ufz,max_u,min_u,max_p,min_p,max_f,min_f
  ELSE
    OPEN(UNIT=441,FILE='fieldInfo.dat',STATUS='REPLACE')
    WRITE(441,20) '   TIME            Qx            Qy            Qz     u_max     u_min     p_max     p_min     f_max     f_min'
    WRITE(441,21) iit,sum_ufx,sum_ufy,sum_ufz,max_u,min_u,max_p,min_p,max_f,min_f
  END IF
  END IF

END SUBROUTINE outputFile

END MODULE io_mod
