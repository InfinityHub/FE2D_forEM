MODULE EM_data_measurements
  ! Author: Jianbo Long, Mar, 2020
  !
  ! Purpose:
  !   set up synthetic measurement positions for modelling
  USE float_precision, ONLY: DPR
  IMPLICIT NONE

  PRIVATE
  ! for synthetic measurements
  INTEGER,PROTECTED  :: nobs = 0
  INTEGER,PROTECTED,ALLOCATABLE  :: matri(:)   ! measurement site attribute
  REAL(DPR),ALLOCATABLE :: xobs(:), yobs(:), zobs(:)

  ! public data
  PUBLIC :: nobs, xobs, yobs, zobs, matri
  ! public procedures
  PUBLIC :: get_synthetic_measurement_positions

CONTAINS

  SUBROUTINE get_synthetic_measurement_positions()
    USE FORTRAN_generic, ONLY: linear_spacing, print_warn_msg, print_general_msg
    USE error_cls
    USE fileio_mod_Peter, ONLY: fileio_exists, fileio_get_fid, fileio_open_read
    USE modelling_parameter, ONLY: modelling

    IMPLICIT NONE
    CHARACTER(LEN=300) :: subroutineTitle = "get_synthetic_measurement_positions"
    INTEGER :: k, null,fid
    TYPE(error_type) :: error

    CALL print_general_msg('getting synthetic measurement sites')
    nobs = 0
    IF(ALLOCATED(xobs))  DEALLOCATE(xobs)
    IF(ALLOCATED(yobs))  DEALLOCATE(yobs)
    IF(ALLOCATED(zobs))  DEALLOCATE(zobs)
    IF(ALLOCATED(matri))  DEALLOCATE(matri)

    ! check if there is a pre-defined measurement site file
    IF( TRIM(ADJUSTL(modelling%measurement_file)) /= '') THEN
       IF( fileio_exists(TRIM(ADJUSTL(modelling%measurement_file))) ) THEN
          ! Open file for reading:
          CALL fileio_get_fid(fid,error)
          IF (error_check(error)) THEN
             CALL error_construct(error,ERROR_GENERAL,'EM_data_measurements',TRIM(ADJUSTL(subroutineTitle)),'failed to get fid ')
             RETURN
          END IF
          CALL fileio_open_read(TRIM(ADJUSTL(modelling%measurement_file)),fid,error)
          IF (error_check(error)) THEN
             CALL error_pass(error,'EM_data_measurements',TRIM(ADJUSTL(subroutineTitle)))
             CALL error_report(error)
             RETURN
          END IF

          ! first line
          READ(fid, *) nobs
          IF( nobs <= 0 ) THEN
             CALL error_construct(error,ERROR_READ,'EM_data_measurements',TRIM(ADJUSTL(subroutineTitle)),'number of measurement sites <= 0!')
             CALL error_report(error)
          END IF
          ALLOCATE( xobs(nobs), yobs(nobs), zobs(nobs), matri(nobs) )

          DO k = 1, nobs
             READ(fid, *) null, xobs(k), zobs(k), matri(k)
          end DO
          yobs = 0.0
          CLOSE(fid)
          RETURN
       end IF
    end IF

    IF( nobs == 0 ) THEN
       CALL print_warn_msg('NO measurement sites are detected for modelling computation !!!')
    END IF
    RETURN
  END SUBROUTINE get_synthetic_measurement_positions
end MODULE EM_Data_Measurements

