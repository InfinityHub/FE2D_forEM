MODULE EMfwd_module
  ! Originally designed for EM problems
  !   Jianbo Long
  IMPLICIT NONE
  PUBLIC :: start_FE_modelling
  PRIVATE
CONTAINS 
  SUBROUTINE start_FE_modelling()
    USE EM_data_measurements, ONLY: get_synthetic_measurement_positions
    USE mesh_concatenate_module, ONLY: mesh_concatenate
    USE num_algorithm_concatenate_differ_operator, ONLY: differ_operator_algorithm_concatenate
    USE linear_system_equation_data, ONLY: assemble_LinearSystem_concatenate,&
         solve_the_linear_system

    USE EM_RHS_concatenate_module, ONLY: EM_RHS_concatenate
    USE EM_data_postprocess, ONLY: data_processing_EM

    USE FORTRAN_generic,ONLY: print_general_msg, print_debug_msg, standard_print_date_and_time
    IMPLICIT NONE

    !CALL standard_print_date_and_time('time information')
    CALL EM_invoke()
    ! check if there are synthetic measurements for post processing (output files)
    CALL get_synthetic_measurement_positions()

    CALL print_general_msg('setting up meshes')
    CALL mesh_concatenate()
    CALL differ_operator_algorithm_concatenate()
    CALL assemble_LinearSystem_concatenate()          
    CALL EM_RHS_concatenate()
    CALL print_general_msg('Solving linear system')
    CALL solve_the_linear_system()
    CALL data_processing_EM()
    RETURN
  end SUBROUTINE start_FE_modelling


  SUBROUTINE EM_invoke()
    ! process inputs for the program, print usage messages.
    USE modelling_parameter, ONLY: reset_output_directory, get_fwd_parameter
    IMPLICIT NONE
    INTEGER :: narg, k, arglen, ierr
    CHARACTER(LEN=256) :: arg, t, inputfile, outDirectory
    !CHARACTER(LEN=*) :: t
    inputfile = ''
    outDirectory = ''
    ! get inputs ( mesh files, model paras, etc) when invoking program
    narg = COMMAND_ARGUMENT_COUNT()
    IF( narg <= 0) THEN
       CALL print_author_info()
       CALL print_usage_info()
       STOP
    end IF

    IF( narg > 0 ) THEN
       DO k = 1, narg
          CALL GET_COMMAND_ARGUMENT(NUMBER=k, VALUE=arg, LENGTH=arglen, STATUS=ierr)  ! get 1st argument
          IF( ierr > 0) THEN
             PRINT*, 'Error: the argument cannot be retrieved: ', TRIM(arg); STOP
          ELSEIF ( ierr == -1 ) THEN
             PRINT*, 'Error: the argument <', TRIM(ADJUSTL(arg)),'> has exceeded the maximum allowance in length (256 characters)!';STOP
          end IF

          t = ADJUSTL(arg)   ! remove leading blanks

          IF( t(1:1) == '-' .AND. arglen > 2 ) THEN
             SELECT CASE( t(2:2) )
             CASE('i')
                ! input file
                inputfile = TRIM( t(3: arglen) )
             CASE('o')
                ! provide a sub-directory for outputs
                outDirectory = TRIM( t(3: arglen) )
             CASE DEFAULT
                PRINT*,'Error: invalid argument flag(s) !'
                PRINT*, TRIM(t)
                CALL print_usage_info();STOP
             END SELECT
          ELSE
             PRINT*,'Error: invalid argument !'
             PRINT*, t
             CALL print_usage_info();STOP
          END IF

       END DO
    END IF

    IF(LEN_TRIM(inputfile) == 0) THEN
       PRINT*,'Error: missing the input file ! Please select an input file.'
       CALL print_usage_info();STOP
    END IF
    CALL get_fwd_parameter(inputfile)
    IF(LEN_TRIM(outDirectory) /= 0) CALL reset_output_directory( outDirectory )

    RETURN
  CONTAINS
    SUBROUTINE print_author_info()
      IMPLICIT NONE
      INTEGER :: timevalues(8)
      CHARACTER(LEN=30) :: var
      CHARACTER(LEN=4) :: year
      CALL DATE_AND_TIME(VALUES=timevalues)
      WRITE(year, '(I4.4)')  timevalues(1)
      PRINT*,'         General 2-D FE modelling program. (unstructured meshes)      '
      PRINT*,'        Jianbo Long, Department of Earth Sciences, MUN. 2017-'//year//'.    '
      PRINT*,'----------------------------------------------------------------------'
      RETURN
    end SUBROUTINE print_author_info

    SUBROUTINE print_usage_info()
      IMPLICIT NONE
      PRINT*,''
      PRINT*,'Usage: FE2D -i<arg> [-o<arg>]'
      PRINT*,''
      PRINT*,'     -i  Select an input file that provides the paths and names of all necessary parameter files, which in turn include '
      PRINT*,'         all set-up parameters (model, numerical method, data output, etc) for the forward modelling. If the file '
      PRINT*,'         is elsewhere, a relative or absolute path can be included.'
      PRINT*,'     -o  Optional. Specify a NAME (not a path) of additional sub-directory for writing all output files from the modelling. '
      PRINT*,'         If not supplied, the default output directory/path, which is set in the parameter files, will be used. If the '
      PRINT*,'         supplied sub-directory does not exist, it will be created.'
      RETURN
    end SUBROUTINE print_usage_info
  end SUBROUTINE EM_invoke  
end MODULE EMfwd_module

PROGRAM FE_FORWARD
  USE EMfwd_module, ONLY: start_FE_modelling
  USE mpi
  IMPLICIT NONE
  INTEGER :: ierr
  CALL MPI_INIT(ierr)
  CALL start_FE_modelling()
  CALL MPI_FINALIZE(ierr)
end PROGRAM FE_FORWARD

