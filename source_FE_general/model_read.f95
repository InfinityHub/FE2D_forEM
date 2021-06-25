MODULE modelling_parameter
  ! read all modelling parameters from files
  USE float_precision, ONLY: DPR
  USE strings_mod, ONLY: WORDLEN
  IMPLICIT NONE
  PRIVATE
  INTEGER,PARAMETER :: PATHLEN = 400
  CHARACTER(LEN=WORDLEN),PROTECTED :: thisModule = "modelling_parameter"
  ! define a type that contains necessary parameters
  TYPE modelling_para_type
     CHARACTER(LEN=2) :: DataType        ! data type, 'CM' (complex) or 'RE' (real)
     ! Output file control
     CHARACTER(LEN=PATHLEN):: measurement_file
     CHARACTER(LEN=PATHLEN):: outputfilepath
     !the reference path/directory is the true parent directory
     CHARACTER(LEN=PATHLEN):: reference_path
     ! Finite element
     LOGICAL :: global_FE_search
     INTEGER :: FEdegree     ! so far: 1 (linear) or 2 (quadratic)
     ! linear system solver
     LOGICAL  :: Iter_solver      !  1: Iterative solvers; otherwise direct solvers
     CHARACTER(LEN=8)   :: iter_solver_name  ! GMRES, BCGSTAB
     LOGICAL  :: solver_verbose              ! if output info
  end type modelling_para_type

  TYPE mesh_para_type
     ! Mesh
     INTEGER :: nregion
     CHARACTER(LEN=PATHLEN) :: regionMarkFile   ! mesh region markers
     INTEGER,ALLOCATABLE :: regionMark(:)
     CHARACTER(LEN=PATHLEN):: meshfilepath
     CHARACTER(LEN=PATHLEN)  :: basefilename   ! base name of mesh files
  end type mesh_para_type

  TYPE :: EM_generic_type
     ! EM general
     REAL(DPR) :: freq
     REAL(DPR) :: omega
     CHARACTER(LEN=PATHLEN)  :: condRegionFile   ! files containing cond list
     INTEGER   ::  ncond
     REAL(DPR),ALLOCATABLE :: cond(:)   ! conductivities of uniform regions
  end type EM_generic_type

  TYPE fwdinp_type
     ! fwd type, pointers to all other input file names.
     CHARACTER(LEN=PATHLEN) :: input_path
     CHARACTER(LEN=PATHLEN) :: input_mesh
     CHARACTER(LEN=PATHLEN) :: input_modelling
     CHARACTER(LEN=PATHLEN) :: input_EM_generic
  end type fwdinp_type

  TYPE(modelling_para_type) :: modelling
  TYPE(mesh_para_type) :: mesh_par
  TYPE(EM_generic_type) :: EM_generic

  PUBLIC :: modelling, mesh_par, EM_generic, get_fwd_parameter, &
       reset_output_directory

CONTAINS

  SUBROUTINE modelling_para_clear(modelling_para)
    IMPLICIT NONE
    TYPE(modelling_para_type),INTENT(INOUT) :: modelling_para
    modelling_para%DataType = ''
    modelling_para%measurement_file = ''
    modelling_para%outputfilepath = ''
    modelling_para%reference_path = ''
    modelling_para%FEdegree = 0
    modelling_para%global_FE_search = .FALSE.
    modelling_para%Iter_solver = .FALSE.
    modelling_para%iter_solver_name = 'NULL'
    modelling_para%solver_verbose = .FALSE.
  end SUBROUTINE modelling_para_clear
  ! ---------------------------------------------
  ! ---------------------------------------------

  SUBROUTINE modelling_para_print(modelling_para)
    IMPLICIT NONE
    TYPE(modelling_para_type),INTENT(IN) :: modelling_para
    PRINT*, "modelling_para%DataType          = ", modelling_para%DataType
    PRINT*, "modelling_para%FEdegree          = ", modelling_para%FEdegree
    PRINT*, "modelling_para%global_FE_search  = ", modelling_para%global_FE_search
    PRINT*, ''
    PRINT*, "modelling_para%Iter_solver       = ", modelling_para%Iter_solver
    IF( modelling_para%Iter_solver ) THEN
       PRINT*, "modelling_para%iter_solver_name  = ", modelling_para%iter_solver_name
    END IF
    PRINT*, "modelling_para%solver_verbose    = ", modelling_para%solver_verbose
  end SUBROUTINE Modelling_Para_Print

  ! ---------------------------------------------
  ! ---------------------------------------------

  SUBROUTINE mesh_para_clear(mesh_para)
    IMPLICIT NONE
    TYPE(mesh_para_type),INTENT(INOUT) :: mesh_para
    mesh_para%nregion = 0
    IF(ALLOCATED(mesh_para%regionMark)) DEALLOCATE(mesh_para%regionMark)
    mesh_para%regionMarkFile = ''
    mesh_para%meshfilepath = ''
    mesh_para%basefilename = ''
  end SUBROUTINE mesh_para_clear

  SUBROUTINE mesh_para_print(mesh_para)
    IMPLICIT NONE
    TYPE(mesh_para_type),INTENT(IN) :: mesh_para
    PRINT*, "mesh_para%nregion        = ", mesh_para%nregion
    PRINT*, "mesh_para%regionMarkFile = ", TRIM(mesh_para%regionMarkFile)
    PRINT*, "mesh_para%regionMark     = ", mesh_para%regionMark
    PRINT*, "mesh_para%meshfilepath   = ", TRIM(mesh_para%meshfilepath)
    PRINT*, "mesh_para%basefilename   = ", TRIM(mesh_para%basefilename)
  end SUBROUTINE mesh_para_print

  ! ---------------------------------------------
  ! ---------------------------------------------

  SUBROUTINE EM_generic_para_clear(EM_para)
    IMPLICIT NONE
    TYPE(EM_generic_type),INTENT(INOUT) :: EM_para
    EM_para%freq = 0.0
    EM_para%omega = 0.0
    EM_para%condRegionFile = ''
    EM_para%ncond = 0
    IF(ALLOCATED(EM_para%cond)) DEALLOCATE(EM_para%cond)
  end SUBROUTINE  EM_generic_para_clear
  ! ---------------------------------------------
  SUBROUTINE EM_generic_para_print(EM_para)
    IMPLICIT NONE
    TYPE(EM_generic_type),INTENT(IN) :: EM_para
    PRINT*, "EM_para%freq           = ", REAL(EM_para%freq)
    PRINT*, "EM_para%condRegionFile = ", TRIM(EM_para%condRegionFile)
    PRINT*, "EM_para%ncond          = ", EM_para%ncond
    PRINT*, "EM_para%cond           = ", REAL(EM_para%cond)
  end SUBROUTINE  EM_generic_para_print
  ! ---------------------------------------------

  SUBROUTINE fwd_clear(fwdinp)
    IMPLICIT NONE
    TYPE(fwdinp_type),INTENT(INOUT) :: fwdinp
    fwdinp%input_path = ''
    fwdinp%input_modelling = ''
    fwdinp%input_mesh = ''
    fwdinp%input_EM_generic = ''
  end SUBROUTINE fwd_clear
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

  SUBROUTINE modelling_read_parameter(modelling,t1,t2,ok,error)
    ! Based on a subroutine from Peter
    ! Deals with a line from an input file in format "name value".
    USE error_cls
    USE strings_mod, ONLY: concatenate
    USE FORTRAN_string_manipulation, ONLY: string_convert_to_lowercase
    USE fileio_mod_peter, ONLY: LINELEN

    IMPLICIT NONE
    TYPE(modelling_para_type), INTENT(INOUT) :: modelling
    CHARACTER(LEN=*), INTENT(INOUT) :: t1,t2 ! hold the name and value strings
    LOGICAL, INTENT(OUT) :: ok ! false if the parameter name is not recognized
    TYPE(error_type), INTENT(INOUT) :: error
    INTEGER :: ierr
    CHARACTER(LEN=LINELEN) :: t

    ! Initialize the ok and ierr variables:
    ok = .TRUE.
    ierr = 0

    ! Any line starting with "#",'!','%' is a commented line that should be skipped:
    t1 = ADJUSTL(t1)
    IF (t1(1:1)=='#') RETURN
    IF (t1(1:1)=='!') RETURN
    IF (t1(1:1)=='%') RETURN

    CALL string_convert_to_lowercase(t1)

    ! Check for parameter name and extract parameter value:
    IF (t1(1:9)=='datatype ') THEN
       modelling%DataType = t2
    ELSE IF (t1(1:17)=='measurement_file ') THEN
       IF (LEN_TRIM(t2)==0) THEN
          CALL error_80()
          IF (error_check(error)) RETURN
       END IF
       modelling%measurement_file = t2       
    ELSE IF (t1(1:15)=='outputfilepath ') THEN
       IF (LEN_TRIM(t2)==0) THEN
          CALL error_80()
          IF (error_check(error)) RETURN
       END IF
       modelling%outputfilepath = t2
       modelling%reference_path = t2

    ELSE IF (t1(1:9)=='fedegree ') THEN
       READ(t2,FMT=*,IOSTAT=ierr) modelling%FEdegree
    ELSE IF (t1(1:25)=='global_fe_search ') THEN
       IF (TRIM(ADJUSTL(t2))== 'f' .OR. TRIM(ADJUSTL(t2))== 'F') THEN
          modelling%global_FE_search = .FALSE.
       ELSEIF (TRIM(ADJUSTL(t2))== 't' .OR. TRIM(ADJUSTL(t2))== 'T') THEN
          modelling%global_FE_search = .TRUE.
       ELSE
          CALL error_70()
          IF (error_check(error)) RETURN
       END IF

    ELSE IF (t1(1:12)=='iter_solver ') THEN
       IF (TRIM(ADJUSTL(t2))== 'f' .OR. TRIM(ADJUSTL(t2))== 'F') THEN
          modelling%Iter_solver = .FALSE.
       ELSEIF (TRIM(ADJUSTL(t2))== 't' .OR. TRIM(ADJUSTL(t2))== 'T') THEN
          modelling%Iter_solver = .TRUE.
       ELSE
          CALL error_70()
          IF (error_check(error)) RETURN
       END IF

    ELSE IF (t1(1:21)=='iter_solver_name ') THEN
       modelling%iter_solver_name = t2
    ELSE IF (t1(1:15)=='solver_verbose ') THEN
       IF (TRIM(ADJUSTL(t2))== 'f' .OR. TRIM(ADJUSTL(t2))== 'F') THEN
          modelling%solver_verbose = .FALSE.
       ELSEIF (TRIM(ADJUSTL(t2))== 't' .OR. TRIM(ADJUSTL(t2))== 'T') THEN
          modelling%solver_verbose = .TRUE.
       ELSE
          CALL error_70()
          IF (error_check(error)) RETURN
       END IF
    ELSE
       ok = .FALSE.
    END IF

    ! Check for error reading parameter value:
    IF (ierr/=0) THEN
       WRITE(t2,FMT='(A,A)') 'parameter ', TRIM(t1)
       CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'modelling_read_parameter',t2)
       !RETURN
    END IF

    RETURN
    ! Error messaging code is below:
  CONTAINS
    SUBROUTINE error_80()
      IMPLICIT NONE
      CALL concatenate(error,t,'empty string for parameter ',TRIM(t1),trimflag_op=.FALSE.)
      IF (error_check(error)) THEN
         CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'modelling_read_parameter')
         RETURN
      END IF
      CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'modelling_read_parameter',t)
    END SUBROUTINE error_80

    SUBROUTINE error_70()
      IMPLICIT NONE
      CALL concatenate(error,t,'unable to recognize the value string for parameter ',TRIM(t1),trimflag_op=.FALSE.)
      IF (error_check(error)) THEN
         CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'modelling_read_parameter')
         RETURN
      END IF
      CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'modelling_read_parameter',t)
    END SUBROUTINE error_70
  END SUBROUTINE modelling_read_parameter

  ! ---------------------------------------------
  ! ---------------------------------------------

  SUBROUTINE mesh_read_parameter(mesh,t1,t2,ok,error)
    ! setting parameters for mesh
    USE error_cls
    USE strings_mod, ONLY: concatenate
    USE FORTRAN_string_manipulation, ONLY: string_convert_to_lowercase
    USE fileio_mod_peter, ONLY: LINELEN

    IMPLICIT NONE
    TYPE(mesh_para_type), INTENT(INOUT) :: mesh
    CHARACTER(LEN=*), INTENT(INOUT) :: t1,t2 ! hold the name and value strings
    LOGICAL, INTENT(OUT) :: ok ! false if the parameter name is not recognized
    TYPE(error_type), INTENT(INOUT) :: error

    INTEGER :: ierr
    CHARACTER(LEN=LINELEN) :: t

    ! Initialize the ok and ierr variables:
    ok = .TRUE.
    ierr = 0

    ! Any line starting with "#",'!','%' is a commented line that should be skipped:
    t1 = ADJUSTL(t1)
    IF (t1(1:1)=='#') RETURN
    IF (t1(1:1)=='!') RETURN
    IF (t1(1:1)=='%') RETURN

    CALL string_convert_to_lowercase(t1)

    ! Check for parameter name and extract parameter value:
    IF (t1(1:8)=='nregion ') THEN
       READ(t2,FMT=*,IOSTAT=ierr) mesh%nregion
    ELSE IF (t1(1:15)=='regionmarkfile ') THEN
       IF (LEN_TRIM(t2)==0) THEN
          CALL error_80()
          IF (error_check(error)) RETURN
       END IF
       mesh%regionMarkFile = t2
    ELSE IF (t1(1:13)=='meshfilepath ') THEN
       IF (LEN_TRIM(t2)==0) THEN
          CALL error_80()
          IF (error_check(error)) RETURN
       END IF
       mesh%meshfilepath = t2
    ELSE IF (t1(1:13)=='basefilename ') THEN
       IF (LEN_TRIM(t2)==0) THEN
          CALL error_80()
          IF (error_check(error)) RETURN
       END IF
       mesh%basefilename = t2
    ELSE
       ok = .FALSE.
    END IF

    ! Check for error reading parameter value:
    IF (ierr/=0) THEN
       WRITE(t2,FMT='(A,A)') 'parameter ', TRIM(t1)
       CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'mesh_read_parameter',t2)
       !RETURN
    END IF

    RETURN
    ! Error messaging code is below:
  CONTAINS
    SUBROUTINE error_80()
      IMPLICIT NONE
      CALL concatenate(error,t,'empty string for parameter ',TRIM(t1),trimflag_op=.FALSE.)
      IF (error_check(error)) THEN
         CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'mesh_read_parameter')
         RETURN
      END IF
      CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'mesh_read_parameter',t)
    END SUBROUTINE error_80
  END SUBROUTINE Mesh_Read_Parameter

  ! ---------------------------------------------
  ! ---------------------------------------------

  SUBROUTINE EM_generic_read_parameter(EM,t1,t2,ok,error)
    ! setting parameters for mesh
    USE error_cls
    USE strings_mod, ONLY: concatenate
    USE FORTRAN_string_manipulation, ONLY: string_convert_to_lowercase
    USE fileio_mod_peter, ONLY: LINELEN

    IMPLICIT NONE
    TYPE(EM_generic_type), INTENT(INOUT) :: EM
    CHARACTER(LEN=*), INTENT(INOUT) :: t1,t2 ! hold the name and value strings
    LOGICAL, INTENT(OUT) :: ok ! false if the parameter name is not recognized
    TYPE(error_type), INTENT(INOUT) :: error

    INTEGER :: ierr
    CHARACTER(LEN=LINELEN) :: t

    ! Initialize the ok and ierr variables:
    ok = .TRUE.
    ierr = 0

    ! Any line starting with "#",'!','%' is a commented line that should be skipped:
    t1 = ADJUSTL(t1)
    IF (t1(1:1)=='#') RETURN
    IF (t1(1:1)=='!') RETURN
    IF (t1(1:1)=='%') RETURN

    CALL string_convert_to_lowercase(t1)

    ! Check for parameter name and extract parameter value:
    IF (t1(1:5)=='freq ') THEN
       READ(t2,FMT=*,IOSTAT=ierr) EM%freq
    ELSE IF (t1(1:15)=='condregionfile ') THEN
       IF (LEN_TRIM(t2)==0) THEN
          CALL error_80()
          IF (error_check(error)) RETURN
       END IF
       EM%condRegionFile = t2
    ELSE
       ok = .FALSE.
    END IF

    ! Check for error reading parameter value:
    IF (ierr/=0) THEN
       WRITE(t2,FMT='(A,A)') 'parameter ', TRIM(t1)
       CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'EM_generic_read_parameter',t2)
       !RETURN
    END IF

    RETURN
    ! Error messaging code is below:
  CONTAINS
    SUBROUTINE error_80()
      IMPLICIT NONE
      CALL concatenate(error,t,'empty string for parameter ',TRIM(t1),trimflag_op=.FALSE.)
      IF (error_check(error)) THEN
         CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'EM_generic_read_parameter')
         RETURN
      END IF
      CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'EM_generic_read_parameter',t)
    END SUBROUTINE error_80
  END SUBROUTINE EM_Generic_Read_Parameter
  ! ---------------------------------------------
  SUBROUTINE fwd_read_parameter(fwdinp,t1,t2,ok,error)
    ! setting parameters (file names) for all input files
    USE error_cls
    USE strings_mod, ONLY: concatenate
    USE FORTRAN_string_manipulation, ONLY: string_convert_to_lowercase
    !USE fileio_mod_peter, ONLY: LINELEN

    IMPLICIT NONE
    TYPE(fwdinp_type), INTENT(INOUT) :: fwdinp
    CHARACTER(LEN=*), INTENT(INOUT) :: t1,t2 ! hold the name and value strings
    LOGICAL, INTENT(OUT) :: ok ! false if the parameter name is not recognized
    TYPE(error_type), INTENT(INOUT) :: error

    INTEGER :: ierr
    !CHARACTER(LEN=LINELEN) :: t

    ! Initialize the ok and ierr variables:
    ok = .TRUE.
    ierr = 0

    ! Any line starting with "#",'!','%' is a commented line that should be skipped:
    t1 = ADJUSTL(t1)
    IF (t1(1:1)=='#') RETURN
    IF (t1(1:1)=='!') RETURN
    IF (t1(1:1)=='%') RETURN

    CALL string_convert_to_lowercase(t1)

    ! Check for parameter name and extract parameter value:
    IF (t1(1:11)=='input_path ') THEN
       fwdinp%input_path = t2
    ELSE IF (t1(1:16)=='input_modelling ') THEN
       fwdinp%input_modelling = t2
    ELSE IF (t1(1:11)=='input_mesh ') THEN
       fwdinp%input_mesh = t2
    ELSE IF (t1(1:17)=='input_em_generic ') THEN
       fwdinp%input_EM_generic = t2
    ELSE
       ok = .FALSE.
    END IF

    ! Check for error reading parameter value:
    IF (ierr/=0) THEN
       WRITE(t2,FMT='(A,A)') 'parameter ', TRIM(t1)
       CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'fwd_read_parameter',t2)
       !RETURN
    END IF
    RETURN
    ! Error messaging code is below:
!!$  CONTAINS
!!$    SUBROUTINE error_80()
!!$      IMPLICIT NONE
!!$      CALL concatenate(error,t,'empty string for parameter ',TRIM(t1),trimflag_op=.FALSE.)
!!$      IF (error_check(error)) THEN
!!$         CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'fwd_read_parameter')
!!$         RETURN
!!$      END IF
!!$      CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'fwd_read_parameter',t)
!!$    END SUBROUTINE error_80
  END SUBROUTINE Fwd_Read_Parameter


  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

  SUBROUTINE modelling_read_inputfile(modelling, inputfile, error)
    USE error_cls
    USE FORTRAN_fileio_mod, ONLY : fileio_exists
    USE fileio_mod_Peter, ONLY: read_input_line, fileio_open_read, fileio_get_fid
    IMPLICIT NONE
    TYPE(modelling_para_type), INTENT(INOUT) :: modelling
    CHARACTER(LEN=*), INTENT(IN) :: inputfile
    TYPE(error_type), INTENT(INOUT) :: error

    LOGICAL :: ok
    INTEGER :: fid,ierr
    CHARACTER(LEN=wordlen) :: t1,t2

    ! Check file exists:
    IF(.NOT.fileio_exists(inputfile)) THEN
       CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'modelling_read_inputfile',&
            'failed to find input file: '//TRIM(inputfile))
       RETURN
    END IF

    ! Clear modelling object:
    CALL modelling_para_clear(modelling)

    ! Open file for reading:
    CALL fileio_get_fid(fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'modelling_read_inputfile')
       RETURN
    END IF
    CALL fileio_open_read(inputfile,fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'modelling_read_inputfile')
       RETURN
    END IF

    ! Loop over each line in the file:
    DO
       ! Read line from file:
       CALL read_input_line(fid,t1,t2,ierr,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'modelling_read_inputfile')
          RETURN
       END IF
       IF (ierr<0) EXIT ! reached end of file
       IF (ierr>0) CYCLE ! commented or empty line
       ! Deal with current line:
       CALL modelling_read_parameter(modelling,t1,t2,ok,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'modelling_read_inputfile')
          RETURN
       END IF
       IF (.NOT.ok) THEN
          WRITE(t2,FMT='(A,A)') 'unknown parameter name: ', TRIM(t1)
          CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'modelling_read_inputfile',t2)
          RETURN
       END IF
    END DO
    ! Close input file:
    CLOSE(UNIT=fid)
    RETURN
  end SUBROUTINE modelling_read_inputfile

  ! ---------------------------------------------
  ! ---------------------------------------------

  SUBROUTINE mesh_read_inputfile(mesh, inputfile, error)
    USE error_cls
    USE fileio_mod_Peter, ONLY: fileio_exists, read_input_line, fileio_open_read, fileio_get_fid
    IMPLICIT NONE
    TYPE(mesh_para_type), INTENT(INOUT) :: mesh
    CHARACTER(LEN=*), INTENT(IN) :: inputfile
    TYPE(error_type), INTENT(INOUT) :: error

    LOGICAL :: ok
    INTEGER :: fid,ierr
    CHARACTER(LEN=wordlen) :: t1,t2

    ! Check file exists:
    IF(.NOT.fileio_exists(inputfile)) THEN
       CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'mesh_read_inputfile','failed to find input file: '//TRIM(ADJUSTL(inputfile)))
       RETURN
    END IF

    ! Clear mesh object:
    CALL mesh_para_clear(mesh)

    ! Open file for reading:
    CALL fileio_get_fid(fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'mesh_read_inputfile')
       RETURN
    END IF
    CALL fileio_open_read(inputfile,fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'mesh_read_inputfile')
       RETURN
    END IF

    ! Loop over each line in the file:
    DO
       ! Read line from file:
       CALL read_input_line(fid,t1,t2,ierr,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'mesh_read_inputfile')
          RETURN
       END IF
       IF (ierr<0) EXIT ! reached end of file
       IF (ierr>0) CYCLE ! commented or empty line
       ! Deal with current line:
       CALL mesh_read_parameter(mesh,t1,t2,ok,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'mesh_read_inputfile')
          RETURN
       END IF
       IF (.NOT.ok) THEN
          WRITE(t2,FMT='(A,A)') 'unknown parameter name: ', TRIM(t1)
          CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'mesh_read_inputfile',t2)
          RETURN
       END IF
    END DO
    ! Close input file:
    CLOSE(UNIT=fid)
    RETURN
  end SUBROUTINE mesh_read_inputfile
  ! ---------------------------------------------
  ! ---------------------------------------------

  SUBROUTINE mesh_read_regionMark(mesh, error)
    ! read values for the component mesh%regionMark
    USE error_cls
    USE fileio_mod_Peter, ONLY: fileio_exists, read_input_line, fileio_open_read, fileio_get_fid
    IMPLICIT NONE
    TYPE(mesh_para_type), INTENT(INOUT) :: mesh
    TYPE(error_type), INTENT(INOUT) :: error
    !LOGICAL :: ok
    INTEGER :: fid,ierr
    CHARACTER(LEN=wordlen) :: t1,t2, t
    INTEGER :: nr

    ! Check file exists:
    ! check mesh parameter regionMarkFile
    IF( LEN_TRIM(mesh%regionMarkFile) /= 0) THEN
       IF(.NOT.fileio_exists(mesh%regionMarkFile)) THEN
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'mesh_read_regionMark',&
               'failed to find input file: '//TRIM(ADJUSTL(mesh%regionMarkFile)))
          RETURN
       END IF
    END IF

    ! Open file for reading:
    CALL fileio_get_fid(fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'mesh_read_regionMark')
       RETURN
    END IF
    CALL fileio_open_read(TRIM(ADJUSTL(mesh%regionMarkFile)),fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'mesh_read_regionMark')
       RETURN
    END IF

    nr = 0
    ! Loop over each line in the file:
    DO
       ! Read line from file:
       CALL read_input_line(fid,t1,t2,ierr,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'mesh_read_regionMark')
          RETURN
       END IF
       IF (ierr<0) EXIT ! reached end of file
       IF (ierr>0) CYCLE ! commented or empty line
       ! Deal with current line:
       IF( TRIM(ADJUSTL(t1)) == 'L') THEN
          READ(t2,FMT=*,IOSTAT=ierr) mesh%nregion
          IF (ierr/=0) THEN
             WRITE(t,FMT='(A,A,A,A)') ':Parameter ', TRIM(t2), ' is not an integer value in the line record: ', TRIM(t1)//'   '//TRIM(t2)
             CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'mesh_read_regionMark',t)
             RETURN
          END IF
          ALLOCATE( mesh%regionMark(mesh%nregion) )
       ELSE
          IF(.NOT. ALLOCATED(mesh%regionMark) ) THEN
             WRITE(t,FMT='(A,A)') 'the line record: ', TRIM(t1)//'   '//TRIM(t2)
             CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'mesh_read_regionMark',t)
             RETURN
          END IF
          nr = nr + 1
          READ(t2,FMT=*,IOSTAT=ierr) mesh%regionMark(nr)
          IF (ierr/=0) THEN
             WRITE(t,FMT='(A,A,A,A)') ':Parameter ', TRIM(t2), ' is not an integer value in the line record: ', TRIM(t1)//'   '//TRIM(t2)
             CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'mesh_read_regionMark',t)
             RETURN
          END IF
       END IF
    END DO
    ! Close input file:
    CLOSE(UNIT=fid)
    RETURN
  end SUBROUTINE Mesh_Read_RegionMark

  ! ---------------------------------------------
  ! ---------------------------------------------

  SUBROUTINE EM_generic_read_inputfile(EM_generic, inputfile, error)
    USE error_cls
    USE fileio_mod_Peter, ONLY: fileio_exists, read_input_line, fileio_open_read, fileio_get_fid
    IMPLICIT NONE
    TYPE(EM_generic_type), INTENT(INOUT) :: EM_generic
    CHARACTER(LEN=*), INTENT(IN) :: inputfile
    TYPE(error_type), INTENT(INOUT) :: error

    LOGICAL :: ok
    INTEGER :: fid,ierr
    CHARACTER(LEN=wordlen) :: t1,t2

    ! Check file exists:
    IF(.NOT.fileio_exists(inputfile)) THEN
       CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'EM_generic_read_inputfile','failed to find input file')
       RETURN
    END IF

    ! Clear modelling object:
    CALL EM_generic_para_clear(EM_generic)

    ! Open file for reading:
    CALL fileio_get_fid(fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'EM_generic_read_inputfile')
       RETURN
    END IF
    CALL fileio_open_read(inputfile,fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'EM_generic_read_inputfile')
       RETURN
    END IF

    ! Loop over each line in the file:
    DO
       ! Read line from file:
       CALL read_input_line(fid,t1,t2,ierr,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'EM_generic_read_inputfile')
          RETURN
       END IF
       IF (ierr<0) EXIT ! reached end of file
       IF (ierr>0) CYCLE ! commented or empty line
       ! Deal with current line:
       CALL EM_generic_read_parameter(EM_generic,t1,t2,ok,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'EM_generic_read_inputfile')
          RETURN
       END IF
       IF (.NOT.ok) THEN
          WRITE(t2,FMT='(A,A)') 'unknown parameter name: ', TRIM(t1)
          CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'EM_generic_read_inputfile',t2)
          RETURN
       END IF
    END DO
    ! Close input file:
    CLOSE(UNIT=fid)
    RETURN
  end SUBROUTINE EM_generic_read_inputfile
  ! ---------------------------------------------------
  ! ---------------------------------------------------
  SUBROUTINE EM_generic_read_conductivity(EM_generic, error)
    ! read values for the component EM_generic%cond
    USE error_cls
    USE fileio_mod_Peter, ONLY: fileio_exists, read_input_line, fileio_open_read, fileio_get_fid
    IMPLICIT NONE
    TYPE(EM_generic_type), INTENT(INOUT) :: EM_generic
    TYPE(error_type), INTENT(INOUT) :: error
    INTEGER :: fid,ierr
    CHARACTER(LEN=wordlen) :: t1,t2, t
    INTEGER :: nr

    ! Check file exists:
    ! check EM_generic parameter condRegionFile
    IF( LEN_TRIM(EM_generic%condRegionFile) /= 0) THEN
       IF(.NOT.fileio_exists(EM_generic%condRegionFile)) THEN
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'EM_generic_read_conductivity',&
               'failed to find input file: '//TRIM(ADJUSTL(EM_generic%condRegionFile)))
          RETURN
       END IF
    END IF

    ! Open file for reading:
    CALL fileio_get_fid(fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'EM_generic_read_conductivity')
       RETURN
    END IF
    CALL fileio_open_read(TRIM(ADJUSTL(EM_generic%condRegionFile)),fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'EM_generic_read_conductivity')
       RETURN
    END IF

    nr = 0
    ! Loop over each line in the file:
    DO
       ! Read line from file:
       CALL read_input_line(fid,t1,t2,ierr,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'EM_generic_read_conductivity')
          RETURN
       END IF
       IF (ierr<0) EXIT ! reached end of file
       IF (ierr>0) CYCLE ! commented or empty line
       ! Deal with current line:
       IF( TRIM(ADJUSTL(t1)) == 'L') THEN
          READ(t2,FMT=*,IOSTAT=ierr) EM_generic%ncond
          IF (ierr/=0) THEN
             WRITE(t,FMT='(A,A,A,A)') ':Parameter ', TRIM(t2), ' is not an integer value in the line record: ', TRIM(t1)//'   '//TRIM(t2)
             CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'EM_generic_read_conductivity',t)
             RETURN
          END IF
          ALLOCATE( EM_generic%cond(EM_generic%ncond) )
       ELSE
          IF(.NOT. ALLOCATED(EM_generic%cond) ) THEN
             WRITE(t,FMT='(A,A)') 'the line record: ', TRIM(t1)//'   '//TRIM(t2)
             CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'EM_generic_read_conductivity',t)
             RETURN
          END IF
          nr = nr + 1
          READ(t2,FMT=*,IOSTAT=ierr) EM_generic%cond(nr)
          IF (ierr/=0) THEN
             WRITE(t,FMT='(A,A,A,A)') ':Parameter ', TRIM(t2), ' is not a numeric value in the line record: ', TRIM(t1)//'   '//TRIM(t2)
             CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'EM_generic_read_conductivity',t)
             RETURN
          END IF
       END IF
    END DO
    ! Close input file:
    CLOSE(UNIT=fid)
    RETURN
  end SUBROUTINE EM_generic_read_conductivity

  SUBROUTINE fwd_read_inputfile(fwdinp, inputfile, error)
    ! read the overall inputfile passed to the modelling program. This file is a 'pointer'.
    USE error_cls
    USE fileio_mod_Peter, ONLY: fileio_exists, read_input_line, fileio_open_read, fileio_get_fid
    IMPLICIT NONE
    TYPE(fwdinp_type), INTENT(INOUT) :: fwdinp
    CHARACTER(LEN=*), INTENT(IN) :: inputfile
    TYPE(error_type), INTENT(INOUT) :: error
    LOGICAL :: ok
    INTEGER :: fid,ierr
    CHARACTER(LEN=wordlen) :: t1,t2

    ! Check file exists:
    IF(.NOT.fileio_exists(inputfile)) THEN
       CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),'fwd_read_inputfile','failed to find input file: <'//TRIM(ADJUSTL(inputfile))//'>')
       RETURN
    END IF

    ! Clear object:
    CALL fwd_clear(fwdinp)

    ! Open file for reading:
    CALL fileio_get_fid(fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'fwd_read_inputfile')
       RETURN
    END IF
    CALL fileio_open_read(inputfile,fid,error)
    IF (error_check(error)) THEN
       CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'fwd_read_inputfile')
       RETURN
    END IF

    ! Loop over each line in the file:
    DO
       ! Read line from file:
       CALL read_input_line(fid,t1,t2,ierr,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'fwd_read_inputfile')
          RETURN
       END IF
       IF (ierr<0) EXIT ! reached end of file
       IF (ierr>0) CYCLE ! commented or empty line
       ! Deal with current line:
       CALL fwd_read_parameter(fwdinp,t1,t2,ok,error)
       IF (error_check(error)) THEN
          CALL error_pass(error,TRIM(ADJUSTL(thisModule)),'fwd_read_inputfile')
          RETURN
       END IF
       IF (.NOT.ok) THEN
          WRITE(t2,FMT='(A,A)') 'unknown parameter name: ', TRIM(t1)
          CALL error_construct(error,ERROR_READ,TRIM(ADJUSTL(thisModule)),'fwd_read_inputfile',t2)
          RETURN
       END IF
    END DO
    ! Close input file:
    CLOSE(UNIT=fid)

    ! check fwd input file path information
    IF( TRIM(ADJUSTL(fwdinp%input_path)) /= '') THEN
       fwdinp%input_modelling = TRIM(ADJUSTL(fwdinp%input_path))//TRIM(ADJUSTL(fwdinp%input_modelling))
       fwdinp%input_mesh = TRIM(ADJUSTL(fwdinp%input_path))//TRIM(ADJUSTL(fwdinp%input_mesh))
       fwdinp%input_EM_generic = TRIM(ADJUSTL(fwdinp%input_path))//TRIM(ADJUSTL(fwdinp%input_EM_generic))
    END IF
    RETURN
  end SUBROUTINE Fwd_Read_Inputfile


  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

  SUBROUTINE set_angular_freq()
    USE constants_module, ONLY: pi
    IMPLICIT NONE
    EM_generic%omega = 2 * pi * EM_generic%freq
    RETURN
  end SUBROUTINE set_angular_freq

  ! ----------------------------------------
  ! ----------------------------------------
  SUBROUTINE reset_output_directory( string_input )
    ! reset the main output destination
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: string_input
    CHARACTER(LEN=500) :: msg
    INTEGER :: ierr
    modelling%reference_path = TRIM(ADJUSTL(modelling%reference_path))//TRIM(ADJUSTL(string_input))//'/'
    ! create the directory on disk if not exist
    CALL EXECUTE_COMMAND_LINE( COMMAND='mkdir -p '//TRIM(ADJUSTL(modelling%reference_path)), CMDSTAT=ierr, CMDMSG=msg )
    IF( ierr /= 0) THEN
       PRINT*, ''; PRINT*,'In <reset_output_directory>'
       PRINT*, "Error in creating the directory <"//TRIM(ADJUSTL(modelling%reference_path))//">. More info : ", TRIM(ADJUSTL(msg))
       STOP
    ENDIF
    RETURN
  end SUBROUTINE reset_output_directory
  ! ----------------------------------------
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
  !-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

  SUBROUTINE get_fwd_parameter(inputfile)
    USE error_cls
    IMPLICIT NONE

    TYPE(fwdinp_type) :: fwdinp
    TYPE(error_type) :: error
    CHARACTER(LEN=*),INTENT(IN) :: inputfile

    CALL modelling_para_clear(modelling)
    CALL mesh_para_clear(mesh_par)
    CALL EM_generic_para_clear(EM_generic)

    ! read the only input file (which points to all other files) for FWD program
    CALL fwd_read_inputfile(fwdinp, inputfile, error)
    IF (error_check(error)) CALL error_report(error)

    ! read modelling parameters
    CALL modelling_read_inputfile(modelling, TRIM(ADJUSTL(fwdinp%input_modelling)), error)
    IF (error_check(error)) CALL error_report(error)

    ! read mesh parameters
    CALL mesh_read_inputfile(mesh_par, TRIM(ADJUSTL(fwdinp%input_mesh)), error)
    IF (error_check(error)) CALL error_report(error)
    CALL mesh_read_regionMark(mesh_par, error)
    IF (error_check(error)) CALL error_report(error)

    ! read EM parameters
    IF(TRIM(ADJUSTL(modelling%DataType)) == 'CM') THEN
       CALL EM_generic_read_inputfile(EM_generic, TRIM(ADJUSTL(fwdinp%input_EM_generic)), error)
       IF (error_check(error)) CALL error_report(error)
       CALL EM_generic_read_conductivity(EM_generic, error)
       IF(error_check(error)) CALL error_report(error)
    END IF

    ! ----- print specs information of modelling to screen ------
    PRINT*,'****************   Finite element modelling  **************'
    PRINT*,''
    PRINT*,'--------  settings of Numerical and physical parameters  ----------------'
    PRINT*, '====> Modelling parameters: '
    CALL modelling_para_print(modelling); PRINT*
    PRINT*, '====> Mesh parameters: '
    CALL mesh_para_print(mesh_par); PRINT*

    IF( TRIM(ADJUSTL(modelling%DataType)) == 'CM' ) THEN
       PRINT*, '====> EM_generic parameters: '
       CALL EM_generic_para_print(EM_generic); PRINT*
    END IF
    PRINT*,'-------------------------------------------------------------------------'
    RETURN
  end SUBROUTINE get_fwd_parameter


end MODULE modelling_parameter
