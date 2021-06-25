MODULE FORTRAN_fileio_mod

  IMPLICIT NONE


CONTAINS

   LOGICAL FUNCTION fileio_exists(filename) RESULT(ok)
   ! Checks if a file exists.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: filename ! name of file to check
      INQUIRE(FILE=TRIM(ADJUSTL(filename)),EXIST=ok) ! PURE functions can't do any I/O
   END FUNCTION fileio_exists  


end MODULE FORTRAN_fileio_mod
