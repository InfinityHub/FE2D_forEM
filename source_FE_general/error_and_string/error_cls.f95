! other
MODULE error_cls
! Class for passing errors up the chain and reporting them.
! Built on a previous version of Peter Lelievre
! Jianbo Long


   IMPLICIT NONE
   
   PRIVATE
   PUBLIC :: error_type, ERROR_GENERAL, ERROR_MEMORY, ERROR_OPENR, ERROR_OPENW, ERROR_READ, ERROR_WRITE
   PUBLIC :: error_construct, error_clear, error_copy, error_check, error_check_and_pass, error_pass, &
        error_report, error_immediate, error_get_desc, error_get_proc, error_debug

   ! Generic procedures:
   INTERFACE error_debug
      MODULE PROCEDURE error_debug_i
      MODULE PROCEDURE error_debug_r
   END INTERFACE error_debug

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! STATIC PROPERTIES
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   INTEGER, PARAMETER :: ERROR_NONE=0, ERROR_INTERNAL=1, ERROR_GENERAL=2, ERROR_MEMORY=4, &
                         ERROR_OPENR=6, ERROR_OPENW=8, ERROR_READ=10, ERROR_WRITE=12
   INTEGER, PARAMETER :: NLEV_MAX=64, LEN_MAX=256
   
!   TYPE e_type
!      INTEGER :: etype
!   END TYPE e_type
!   TYPE(e_type), PARAMETER :: ERROR_NONE   =e_type(0) , ERROR_INTERNAL=e_type(1) , &
!                              ERROR_GENERAL=e_type(2) , ERROR_MEMORY  =e_type(4) , &
!                              ERROR_OPENR  =e_type(6) , ERROR_OPENW   =e_type(8) , &
!                              ERROR_READ   =e_type(10), ERROR_WRITE   =e_type(12)

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! TYPE DEFINITION
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   TYPE error_type
      PRIVATE
      INTEGER :: etype=ERROR_NONE ! type of error (ONE of those defined above)
      INTEGER :: ierr=0 ! e.g. result of a read or allocation statement
      CHARACTER(LEN=LEN_MAX) :: desc='' ! description (may be a file or variable name or a general error message)
      INTEGER :: nlev=0 ! number of levels up which the error message has been passed
      CHARACTER(LEN=132), ALLOCATABLE, DIMENSION(:) :: modules, procedures ! record of where the error has been passed through
   END TYPE error_type

CONTAINS

!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! CONSTRUCT, CLEAR, DEEP COPY
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   PURE SUBROUTINE error_construct(error,etype,modu,proc,desc,ierr_op)
   ! Constructs the error object.
      IMPLICIT NONE
      TYPE(error_type), INTENT(INOUT) :: error
      INTEGER, INTENT(IN) :: etype ! type of error (ONE of the integer values defined in this class)
      CHARACTER(LEN=*), INTENT(IN) :: modu, proc ! module and procedure name where the error occurred
      CHARACTER(LEN=*), INTENT(IN) :: desc ! description
      INTEGER, INTENT(IN), OPTIONAL :: ierr_op ! e.g. result of a read or allocation statement
      INTEGER :: ierr
      ! Check for non-empty error object:
      IF (error%etype/=ERROR_NONE) THEN
         error%etype = ERROR_INTERNAL
         error%desc = 'ERROR: NON-EMPTY ERROR OBJECT SUPPLIED TO ERROR_CONSTRUCT.'
         error%modules(1) = modu
         error%procedures(1) = proc
         RETURN
      END IF
      ! Check for no error etype:
      IF (etype==ERROR_NONE) THEN ! RETURN
         CALL error_clear(error)
         error%etype = ERROR_INTERNAL
         error%desc = 'ERROR: ERROR_NONE TYPE SUPPLIED TO ERROR_CONSTRUCT.'
         RETURN
      END IF
      ! Check for internal error etype:
      IF (etype==ERROR_INTERNAL) THEN ! RETURN
         CALL error_clear(error)
         error%etype = ERROR_INTERNAL
         error%desc = 'ERROR: ERROR_INTERNAL TYPE SUPPLIED TO ERROR_CONSTRUCT.'
         RETURN
      END IF
      ! Allocate:
      ALLOCATE(error%modules(NLEV_MAX),error%procedures(NLEV_MAX),STAT=ierr)
      ! Check for allocation error:
      IF (ierr/=0) THEN
         CALL error_clear(error)
         error%etype = ERROR_INTERNAL
         error%desc = 'ERROR: ALLOCATION FAILED IN ERROR_CONSTRUCT.'
         RETURN
      END IF
      ! Construct:
      error%etype = etype
      IF (PRESENT(ierr_op)) error%ierr = ierr_op
      error%desc = desc
      error%nlev = 1
      error%modules(1) = modu
      error%procedures(1) = proc
   END SUBROUTINE error_construct
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   PURE SUBROUTINE error_clear(error)
   ! Clears the error object.
      IMPLICIT NONE
      TYPE(error_type), INTENT(INOUT) :: error
      error%etype = ERROR_NONE
      error%ierr = 0
      error%desc = ''
      error%nlev = 0
      IF (ALLOCATED(error%modules)) DEALLOCATE(error%modules)
      IF (ALLOCATED(error%procedures)) DEALLOCATE(error%procedures)
   END SUBROUTINE error_clear
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   PURE SUBROUTINE error_copy(error1,error2)
   ! Deep copy for error objects.
      IMPLICIT NONE
      TYPE(error_type), INTENT(IN) :: error1 ! information is copied FROM this object
      TYPE(error_type), INTENT(OUT) :: error2 ! information is copied TO this object (or a new error if an error occurs here)
      INTEGER :: i,ierr
      ! Clear:
      CALL error_clear(error2)
      ! Allocate:
      ALLOCATE(error2%modules(NLEV_MAX),error2%procedures(NLEV_MAX),STAT=ierr)
      ! Check for allocation error:
      IF (ierr/=0) THEN
         CALL error_clear(error2)
         CALL error_construct(error2,ERROR_MEMORY,'error_cls','error_copy','error2 arrays',ierr)
         RETURN
      END IF
      ! Copy:
      error2%etype = error1%etype
      error2%ierr  = error1%ierr
      error2%desc  = error1%desc
      error2%nlev  = error1%nlev
      DO i=1,error1%nlev
         error2%modules(i)    = error1%modules(i)
         error2%procedures(i) = error1%procedures(i)
      END DO
   END SUBROUTINE error_copy
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! GETTERS
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   LOGICAL PURE FUNCTION error_check(error) RESULT(check)
   ! Returns true if the error object holds an error.
      IMPLICIT NONE
      TYPE(error_type), INTENT(IN) :: error
      check = ( error%etype /= ERROR_NONE )
   END FUNCTION error_check

   LOGICAL FUNCTION error_check_and_pass(error,modu,proc) RESULT(check)
   ! combines error_check and error_pass into one line (JIANBO)
      IMPLICIT NONE
      TYPE(error_type) :: error
      CHARACTER(LEN=*), INTENT(IN) :: modu, proc ! module and procedure name where the error occurred
      IF(error_check(error)) THEN
         check = .TRUE.
         CALL error_pass(error, modu, proc)
      ELSE
         check = .FALSE.
      END IF
   END FUNCTION error_check_and_pass   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   PURE SUBROUTINE error_get_desc(error,desc)
   ! Returns the error description.
      IMPLICIT NONE
      TYPE(error_type), INTENT(IN) :: error
      CHARACTER(LEN=*), INTENT(OUT) :: desc
      IF (error%etype==ERROR_NONE) THEN
         desc = ''
      ELSE
         desc = error%desc
      END IF
   END SUBROUTINE error_get_desc
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   PURE SUBROUTINE error_get_proc(error,proc)
   ! Returns the highest error procedure.
      IMPLICIT NONE
      TYPE(error_type), INTENT(IN) :: error
      CHARACTER(LEN=*), INTENT(OUT) :: proc
      IF (error%etype==ERROR_NONE) THEN
         proc = ''
      ELSE
         proc = error%procedures(1)
      END IF
   END SUBROUTINE error_get_proc
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! PUBLIC METHODS
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE error_pass(error,modu,proc,s_op)
   ! Adds new module and procedure names to the pass record.
      IMPLICIT NONE
      TYPE(error_type), INTENT(INOUT) :: error
      CHARACTER(LEN=*), INTENT(IN) :: modu, proc ! module and procedure name where the error occurred
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: s_op ! tacked onto the previous subroutine in the chain
      ! Check for empty error:
      IF (.NOT.error_check(error)) THEN ! RETURN
         CALL error_construct(error,ERROR_GENERAL,modu,proc,'ERROR: EMPTY ERROR OBJECT SUPPLIED TO ERROR_PASS.')
         error%etype = ERROR_INTERNAL
         RETURN
      END IF
      ! Check for internal error:
      IF (error%etype==ERROR_INTERNAL) RETURN ! don't pass in this situation
      ! Check for optional input:
      IF (PRESENT(s_op).AND.(error%nlev>0)) THEN
         ! Tack string onto the previous subroutine in the chain:
         error%procedures(error%nlev) = TRIM(error%procedures(error%nlev)) // TRIM(s_op)
      END IF
      ! Increment the pass record:
      error%nlev = error%nlev + 1
      ! Check for overrunning storage:
      IF (error%nlev>(NLEV_MAX+1)) RETURN ! so that the code in the if block below is only processed once
      IF (error%nlev>NLEV_MAX) THEN ! overwrite first module and procedure record with internal message
         error%modules(1) = 'HIT NLEV_MAX IN ERROR_PASS'
         error%procedures(1) = '(procedure chain has been truncated)'
         RETURN
      END IF
      ! Add passed module and procedure information:
      error%modules(error%nlev) = TRIM(modu)
      error%procedures(error%nlev) = TRIM(proc)
   END SUBROUTINE error_pass
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE error_debug_i(error,p,m,s,i)
   ! Convenience method for debugging.
      IMPLICIT NONE
      TYPE(error_type), INTENT(INOUT) :: error
      CHARACTER(LEN=*), INTENT(IN) :: p,m,s
      INTEGER, DIMENSION(:), INTENT(IN) :: i
      INTEGER :: ierr
      CHARACTER(LEN=LEN_MAX) :: t
      WRITE(UNIT=t,FMT=*,IOSTAT=ierr) TRIM(s),i
      IF (ierr/=0) t(LEN_MAX:LEN_MAX) = '&'
      CALL error_construct(error,ERROR_GENERAL,p,m,t)
   END SUBROUTINE error_debug_i

   PURE SUBROUTINE error_debug_r(error,p,m,s,r)
   ! Convenience method for debugging.
      IMPLICIT NONE
      TYPE(error_type), INTENT(INOUT) :: error
      CHARACTER(LEN=*), INTENT(IN) :: p,m,s
      REAL, DIMENSION(:), INTENT(IN) :: r
      INTEGER :: ierr
      CHARACTER(LEN=LEN_MAX) :: t
      WRITE(UNIT=t,FMT=*,IOSTAT=ierr) TRIM(s),r
      IF (ierr/=0) t(LEN_MAX:LEN_MAX) = '&'
      CALL error_construct(error,ERROR_GENERAL,p,m,t)
   END SUBROUTINE error_debug_r
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! I/O
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   SUBROUTINE error_report(error,fid_op)
   ! Reports the error and kills execution unless fid_op is present.
      IMPLICIT NONE
      TYPE(error_type), INTENT(IN) :: error
      INTEGER, INTENT(IN), OPTIONAL :: fid_op ! optional file id (file should be open already)
      PRINT *, '---- EXECUTION STOPPED DUE TO ERROR (DETAILS BELOW) ----'
      ! Print the message and chain:
      IF (PRESENT(fid_op)) THEN
         CALL print_message(error,fid_op=fid_op)
         CALL print_chain(error,fid_op=fid_op)
         ! Don't kill execution if a file ID is supplied!
      ELSE
         CALL print_message(error)
         CALL print_chain(error)
         ! Kill execution:
         PRINT *, '---- EXECUTION STOPPED DUE TO ERROR (DETAILS ABOVE) ----'
         STOP
      END IF
   END SUBROUTINE error_report
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   SUBROUTINE print_message(error,fid_op)
   ! Prints the error message.
      IMPLICIT NONE
      TYPE(error_type), INTENT(IN) :: error
      INTEGER, INTENT(IN), OPTIONAL :: fid_op
      INTEGER :: fid
      CHARACTER(LEN=LEN_MAX) :: errstr
      ! Check optional input:
      IF (PRESENT(fid_op)) THEN
         fid = fid_op
      ELSE
         fid = -1
      END IF
      ! Check type of error and display message:
      IF (error%nlev<=0) THEN
         IF (fid<=0) THEN
            PRINT *, 'ERROR: UNSPECIFIED ERROR ENCOUNTERED'
         ELSE
            WRITE(UNIT=fid,FMT='(A)') 'ERROR: UNSPECIFIED ERROR ENCOUNTERED'
         END IF
         RETURN
      END IF
      SELECT CASE(error%etype)
         CASE(ERROR_GENERAL)
            CALL generr(error%modules(1),error%procedures(1),error%desc,errstr)
         CASE(ERROR_INTERNAL)
            CALL interr(error%modules(1),error%procedures(1),error%desc,errstr)
         CASE(ERROR_MEMORY)
            CALL memerr(error%modules(1),error%procedures(1),error%desc,errstr)
         CASE(ERROR_OPENR)
            CALL openrerr(error%modules(1),error%procedures(1),error%desc,errstr)
         CASE(ERROR_OPENW)
            CALL openwerr(error%modules(1),error%procedures(1),error%desc,errstr)
         CASE(ERROR_READ)
            CALL readerr(error%modules(1),error%procedures(1),error%desc,errstr)
         CASE(ERROR_WRITE)
            CALL writerr(error%modules(1),error%procedures(1),error%desc,errstr)
         CASE(ERROR_NONE)
            CALL noerr(error%modules(1),error%procedures(1),error%desc,errstr)
         CASE DEFAULT ! ERROR_NONE
            IF (fid<=0) THEN
               PRINT *, 'ERROR: UNKNOWN ERROR TYPE ENCOUNTERED - ',error%etype
               PRINT *, TRIM(error%desc)
            ELSE
               WRITE(UNIT=fid,FMT='(A,1I3)') 'ERROR: UNKNOWN ERROR TYPE ENCOUNTERED - ',error%etype
               WRITE(UNIT=fid,FMT='(A)') TRIM(error%desc)
            END IF
            RETURN
      END SELECT
      IF (fid<=0) THEN
         PRINT *, TRIM(errstr)
      ELSE
         WRITE(UNIT=fid,FMT='(A)') TRIM(errstr)
      END IF
   END SUBROUTINE print_message
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   SUBROUTINE print_chain(error,fid_op)
   ! Prints the error chain.
      IMPLICIT NONE
      TYPE(error_type), INTENT(IN) :: error
      INTEGER, INTENT(IN), OPTIONAL :: fid_op
      INTEGER :: i,fid
      ! Check optional input:
      IF (PRESENT(fid_op)) THEN
         fid = fid_op
      ELSE
         fid = -1
      END IF
      ! Check:
      IF (error%nlev<=1) RETURN ! don't bother printing if there is only ONE level (or none)
      IF ((.NOT.ALLOCATED(error%modules)).OR.(.NOT.ALLOCATED(error%procedures))) THEN
         IF (fid<=0) THEN
            PRINT *, 'ERROR: CHAIN ARRAYS NOT ALLOCATED.'
            PRINT *, 'Calling etype:       ',error%etype
            PRINT *, 'Calling description: ',TRIM(error%desc)
         ELSE
            WRITE(UNIT=fid,FMT='(A)') 'ERROR: CHAIN ARRAYS NOT ALLOCATED.'
            WRITE(UNIT=fid,FMT='(A,1I3)') 'Calling etype:       ',error%etype
            WRITE(UNIT=fid,FMT='(A,A)') 'Calling description: ',TRIM(error%desc)
         END IF
         RETURN
      END IF
      ! Print:
      IF (fid<=0) THEN
         PRINT *, '---- ERROR CHAIN ----'
      ELSE
         WRITE(UNIT=fid,FMT='(A)') '---- ERROR CHAIN ----'
      END IF
      DO i=1,error%nlev
         IF (fid<=0) THEN
            PRINT *, TRIM(error%modules(i))," / ",TRIM(error%procedures(i))
         ELSE
            WRITE(UNIT=fid,FMT='(A,A,A)') TRIM(error%modules(i))," / ",TRIM(error%procedures(i))
         END IF
      END DO
   END SUBROUTINE print_chain
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

   PURE SUBROUTINE noerr(modname,subname,desc,errstr)
   ! "NULL TYPE ERROR in modname/subname: desc"
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: modname ! name of module containing calling subroutine
      CHARACTER(LEN=*), INTENT(IN) :: subname ! name of calling subroutine (i.e. where error occurred)
      CHARACTER(LEN=*), INTENT(IN) :: desc ! description of error
      CHARACTER(LEN=*), INTENT(OUT) :: errstr ! the string to print is dumped into this variable
      errstr = "NULL TYPE ERROR in " // TRIM(modname) // "/" // TRIM(subname) // ": " // TRIM(desc)
   END SUBROUTINE noerr

   PURE SUBROUTINE generr(modname,subname,desc,errstr)
   ! "Error in modname/subname: desc"
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: modname ! name of module containing calling subroutine
      CHARACTER(LEN=*), INTENT(IN) :: subname ! name of calling subroutine (i.e. where error occurred)
      CHARACTER(LEN=*), INTENT(IN) :: desc ! description of error
      CHARACTER(LEN=*), INTENT(OUT) :: errstr ! the string to print is dumped into this variable
      errstr = "Error in " // TRIM(modname) // "/" // TRIM(subname) // ": " // TRIM(desc)
   END SUBROUTINE generr

   PURE SUBROUTINE interr(modname,subname,desc,errstr)
   ! "Internal error in modname/subname: desc"
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: modname ! name of module containing calling subroutine
      CHARACTER(LEN=*), INTENT(IN) :: subname ! name of calling subroutine (i.e. where error occurred)
      CHARACTER(LEN=*), INTENT(IN) :: desc ! description of error
      CHARACTER(LEN=*), INTENT(OUT) :: errstr ! the string to print is dumped into this variable
      errstr = "Error in " // TRIM(modname) // "/" // TRIM(subname) // ": " // TRIM(desc)
   END SUBROUTINE interr

   PURE SUBROUTINE memerr(modname,subname,varname,errstr)
   ! "Memory error in modname/subname while allocating varname"
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: modname ! name of module containing calling subroutine
      CHARACTER(LEN=*), INTENT(IN) :: subname ! name of calling subroutine (i.e. where error occurred)
      CHARACTER(LEN=*), INTENT(IN) :: varname ! name or description of variable being allocated
      CHARACTER(LEN=*), INTENT(OUT) :: errstr ! the string to print is dumped into this variable
      errstr = "Memory error in " // TRIM(modname) // "/" // TRIM(subname) // " while allocating " // TRIM(varname)
   END SUBROUTINE memerr

   PURE SUBROUTINE openrerr(modname,subname,filename,errstr)
   ! "Error in modname/subname while opening filename for reading"
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: modname ! name of module containing calling subroutine
      CHARACTER(LEN=*), INTENT(IN) :: subname ! name of calling subroutine (i.e. where error occurred)
      CHARACTER(LEN=*), INTENT(IN) :: filename ! name or description of file being opened
      CHARACTER(LEN=*), INTENT(OUT) :: errstr ! the string to print is dumped into this variable
      errstr = "Error in " // TRIM(modname) // "/" // TRIM(subname) // " while opening " // TRIM(filename) // " for reading"
   END SUBROUTINE openrerr

   PURE SUBROUTINE openwerr(modname,subname,filename,errstr)
   ! "Error in modname/subname while opening filename for writing"
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: modname ! name of module containing calling subroutine
      CHARACTER(LEN=*), INTENT(IN) :: subname ! name of calling subroutine (i.e. where error occurred)
      CHARACTER(LEN=*), INTENT(IN) :: filename ! name or description of file being opened
      CHARACTER(LEN=*), INTENT(OUT) :: errstr ! the string to print is dumped into this variable
      errstr = "Error in " // TRIM(modname) // "/" // TRIM(subname) // " while opening " // TRIM(filename) // " for writing"
   END SUBROUTINE openwerr
   
   PURE SUBROUTINE readerr(modname,subname,varname,errstr)
   ! "Error in modname/subname while reading varname"
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: modname ! name of module containing calling subroutine
      CHARACTER(LEN=*), INTENT(IN) :: subname ! name of calling subroutine (i.e. where error occurred)
      CHARACTER(LEN=*), INTENT(IN) :: varname ! name or description of item being read
      CHARACTER(LEN=*), INTENT(OUT) :: errstr ! the string to print is dumped into this variable
      errstr = "Error in " // TRIM(modname) // "/" // TRIM(subname) // " while reading " // TRIM(varname)
   END SUBROUTINE readerr

   PURE SUBROUTINE writerr(modname,subname,varname,errstr)
   ! "Error in modname/subname while writing varname"
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: modname ! name of module containing calling subroutine
      CHARACTER(LEN=*), INTENT(IN) :: subname ! name of calling subroutine (i.e. where error occurred)
      CHARACTER(LEN=*), INTENT(IN) :: varname ! name or description of item being written
      CHARACTER(LEN=*), INTENT(OUT) :: errstr ! the string to print is dumped into this variable
      errstr = "Error in " // TRIM(modname) // "/" // TRIM(subname) // " while writing " // TRIM(varname)
   END SUBROUTINE writerr
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
! STATIC METHODS
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
   
   SUBROUTINE error_immediate(etype,modu,proc,desc,ierr_op)
   ! Immediately reports an error and kills execution.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: etype ! type or error (ONE of the integer values defined in this class)
      CHARACTER(LEN=*), INTENT(IN) :: modu, proc ! module and procedure name where the error occurred
      CHARACTER(LEN=*), INTENT(IN) :: desc ! description
      INTEGER, INTENT(IN), OPTIONAL :: ierr_op ! e.g. result of a read or allocation statement
      TYPE(error_type) :: error
      ! Construct error object:
      IF (PRESENT(ierr_op)) THEN
         CALL error_construct(error,etype,modu,proc,desc,ierr_op)
      ELSE
         CALL error_construct(error,etype,modu,proc,desc)
      END IF
      ! Report:
      CALL error_report(error)
   END SUBROUTINE error_immediate
   
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!
!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!-----!

END MODULE error_cls
