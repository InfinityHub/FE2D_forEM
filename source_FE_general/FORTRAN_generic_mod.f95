MODULE FORTRAN_generic
  ! SOME general-purpose subroutines related to array manipulation


  ! Jianbo Long, April, 2018
  USE float_precision, ONLY: DPR, CDPR

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: linear_spacing, natural_log_spacing, meshing_grid_2D

  PUBLIC :: add_to_list_real, add_to_list_integer, sort_array_integer,sort_array_real,&
       whether_two_integer_arrays_identical, complex_dot_product

  PUBLIC :: print_general_msg, print_warn_msg, print_error_msg, print_debug_msg, &
       standard_print_date_and_time, percentage_display, find_intersect_integer


  INTERFACE print_debug_msg
     MODULE PROCEDURE print_debug_msg_r, print_debug_msg_i, print_debug_msg_c
  END INTERFACE

CONTAINS

  ! ==========================================================================================

  subroutine add_to_list_real(TargetList, element)
    !! To grow the size of an array on-the-fly; each time add one new entry
    !!    to the previous list/array.
    implicit none
    real(DPR),allocatable:: TargetList(:), newlist(:)
    real(DPR)     element
    INTEGER     isize

    if( allocated(TargetList) ) then
       isize = size(TargetList)

       allocate(newlist(isize + 1) )

       newlist(1: isize) = TargetList
       newlist(isize + 1) = element
       deallocate( TargetList)
       call move_alloc(newlist, TargetList)
    else
       allocate( TargetList(1) )
       TargetList(1) = element
    end if

    return
  end subroutine add_to_list_real

  ! ==========================================================================================

  subroutine add_to_list_integer(TargetList, element)
    !! To grow the size of an array on-the-fly; each time add one new entry
    !!    to the previous list/array.
    implicit none
    integer,allocatable:: TargetList(:), newlist(:)
    integer     element, isize

    if( allocated(TargetList) ) then
       isize = size(TargetList)
       allocate(newlist(isize + 1) )

       newlist(1: isize) = TargetList
       newlist(isize + 1) = element
       deallocate( TargetList)
       call move_alloc(newlist, TargetList)
    else
       allocate( TargetList(1) )
       TargetList(1) = element
    end if

  end subroutine add_to_list_integer


  ! ==========================================================================================

  FUNCTION sort_array_integer(array) RESULT(newlist)
    !! PROCEDURE sort_array_integer() sorts integer array elements in the ascending order. It uses
    !   intrinsic function maxloc/minloc, maxval/minval.
    !
    ! -- Note: there is no need to assume that all elements are distinct !
    !
    ! -- minloc/maxloc returns an rank-one array with size of one, i.e., a(1), if the input is a rank-1 array
    IMPLICIT NONE
    !INTEGER,INTENT(IN),ALLOCATABLE  :: array(:)
    INTEGER,INTENT(IN) :: array(:)
    ! INTEGER  :: newlist( SIZE(array) )
    INTEGER,ALLOCATABLE  :: newlist(: )

    INTEGER,ALLOCATABLE :: temparray(:), shrink(:)
    INTEGER :: i, k,m, isize



    isize = SIZE(array)
    IF( isize < 1) CALL report_error('--Inputed array has less than ONE element--')

    IF( ALLOCATED(newlist) )  DEALLOCATE(newlist); ALLOCATE( newlist(isize) )

    IF( isize == 1 )  THEN
       newlist = array
    ELSE

       IF( ALLOCATED(temparray) )  DEALLOCATE(temparray); ALLOCATE( temparray(isize) )
       temparray = array

       DO i = 1, isize - 1

          newlist(i) = temparray( MINLOC(temparray, DIM = 1) )  !! temparray( MINLOC(temparray) )

          IF( SIZE(temparray) == 2 )  THEN   !! when i == isize - 1
             newlist( isize ) = MAXVAL(temparray)
             EXIT
          ENDIF

          IF( ALLOCATED(shrink) )  DEALLOCATE(shrink);  ALLOCATE( shrink( isize - i ) )

          !! shrink = PACK(temparray, temparray /= newlist(i) )
          m = 0
          DO k = 1, SIZE(temparray)
             IF( k /= MINLOC(temparray, DIM = 1) ) THEN
                m = m + 1;  shrink(m) = temparray(k)
             ENDIF
          ENDDO

          IF( ALLOCATED(temparray) )  DEALLOCATE(temparray);  ALLOCATE( temparray( isize - i ) )
          temparray = shrink

       END DO

    END IF

    RETURN
  end FUNCTION sort_array_integer


  ! ==========================================================================================

  FUNCTION sort_array_real(array) RESULT(newlist)
    !! PROCEDURE sort_array_real() sorts real array elements in the ascending order. It uses
    !   intrinsic function maxloc/minloc, maxval/minval.
    !
    ! -- Note: there is no need to assume that all elements are distinct !
    !
    ! -- minloc/maxloc returns an rank-one array with size of one, i.e., a(1), if the input is a rank-1 array
    IMPLICIT NONE
    REAL(DPR),INTENT(IN),ALLOCATABLE  :: array(:)

    REAL(DPR),ALLOCATABLE  :: newlist(:)

    REAL(DPR),ALLOCATABLE :: temparray(:), shrink(:)
    INTEGER :: i, k,m, isize



    isize = SIZE(array)
    IF( isize < 1) CALL report_error('--Inputed array has less than ONE element--')

    IF( ALLOCATED(newlist) )  DEALLOCATE(newlist); ALLOCATE( newlist(isize) )

    IF( isize == 1 )  THEN
       newlist = array
    ELSE

       IF( ALLOCATED(temparray) )  DEALLOCATE(temparray); ALLOCATE( temparray(isize) )
       temparray = array

       DO i = 1, isize - 1

          newlist(i) = temparray( MINLOC(temparray, DIM = 1) )  !! temparray( MINLOC(temparray) )

          IF( SIZE(temparray) == 2 )  THEN   !! when i == isize - 1
             newlist( isize ) = MAXVAL(temparray)
             EXIT
          ENDIF

          IF( ALLOCATED(shrink) )  DEALLOCATE(shrink);  ALLOCATE( shrink( isize - i ) )

          !! shrink = PACK(temparray, temparray /= newlist(i) )
          m = 0
          DO k = 1, SIZE(temparray)
             IF( k /= MINLOC(temparray, DIM = 1) ) THEN
                m = m + 1;  shrink(m) = temparray(k)
             ENDIF
          ENDDO

          IF( ALLOCATED(temparray) )  DEALLOCATE(temparray);  ALLOCATE( temparray( isize - i ) )
          temparray = shrink

       END DO

    END IF

    RETURN
  end FUNCTION sort_array_real

  ! ==========================================================================================
  SUBROUTINE find_intersect_integer(array1, array2, results, errFlag)
    ! find the intersect elements of two (unsorted) integer arrays. Often, the size of the result
    !  is not known before searching.
    !
    ! Note: the size of the returned result can be 0 (no intersect, errFlag will be 1)
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: array1(:), array2(:)
    INTEGER,INTENT(INOUT) :: errFlag
    INTEGER,INTENT(INOUT),ALLOCATABLE :: results(:)

    INTEGER :: k, nin
    LOGICAL,ALLOCATABLE :: keep(:)

    ALLOCATE( keep(SIZE(array1)) );  keep  = .FALSE.

    DO k = 1, SIZE( array1 )
       IF( ANY(array1(k) == array2) )  keep(k) = .TRUE.
    end DO

    IF( ALLOCATED(results) )  DEALLOCATE(results);

    nin = COUNT(keep .EQV. .TRUE.)
    IF( nin > 0 ) THEN
       ALLOCATE( results( nin ) )
       results = PACK( array1, keep)
       errFlag = 0
    ELSE
       ! the case of no intersects
       errFlag = 1

       RETURN
    END IF


    RETURN
  end SUBROUTINE find_intersect_integer




  ! ==========================================================================================
  FUNCTION whether_two_integer_arrays_identical(array1, array2) RESULT(res)
    ! Given two integer arrays with the same dimension, determine if the two arrays contain exactly
    !    the same integers. The two arrays are not required to be sorted upon input.
    !

    IMPLICIT NONE

    INTEGER,INTENT(IN) :: array1(:), array2(:)
    !INTEGER,INTENT(INOUT) :: errFlag
    LOGICAL :: res


!!$    IF( SIZE(array1) == SIZE(array2) )  THEN
!!$       errFlag = 0
!!$    ELSE
!!$       errFlag = 1
!!$    end IF


    IF( ANY( (sort_array_integer(array1) - sort_array_integer(array2)) /= 0 ) )  THEN
       res = .FALSE.
    ELSE
       res = .TRUE.
    END IF


    RETURN
  end FUNCTION whether_two_integer_arrays_identical


  ! ==========================================================================================

  SUBROUTINE standard_print_date_and_time(msg)
    IMPLICIT NONE

    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: msg   ! a printed message
    INTEGER :: timevalues(8)
    CHARACTER(LEN=30) :: var
    CALL DATE_AND_TIME(VALUES=timevalues)

    WRITE(var, '(A,I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I3.3)') ' ',timevalues(1),'/',timevalues(2),&
         '/',timevalues(3), '; ', timevalues(5),':',timevalues(6),':',timevalues(7),':',timevalues(8)
    IF(PRESENT(msg)) CALL print_debug_msg( msg, var )
    RETURN
  end SUBROUTINE standard_print_date_and_time
  ! ---------------------------------

  SUBROUTINE print_general_msg( msg )
    ! write a general message to screen during running time
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: msg
    PRINT*,"--"//TRIM(msg)
    RETURN
  end SUBROUTINE print_general_msg
  ! ---------------------------------
  SUBROUTINE print_warn_msg( msg )
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: msg
    PRINT*,"***"//TRIM(msg)//"***"
    RETURN
  end SUBROUTINE
  ! ---------------------------------
  SUBROUTINE print_error_msg( msg )
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: msg
    PRINT*,"================================================="
    PRINT*,"#### Error:"//TRIM(msg)//"."
    RETURN
  end SUBROUTINE
  ! ---------------------------------
  SUBROUTINE print_debug_msg_r( msg, a )
    ! write a dubug message to screen during running time
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: msg
    REAL(DPR),INTENT(IN) :: a
    !CHARACTER(LEN=256) :: var
    !WRITE(var,FMT="(A)") a
    PRINT*,"---DEBUG---"//TRIM(msg)//"=", a 
    RETURN
  end SUBROUTINE print_debug_msg_r
  ! ---------------------------------
  SUBROUTINE print_debug_msg_i( msg, a )
    ! write a dubug message to screen during running time
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: msg
    INTEGER,INTENT(IN) :: a
    !CHARACTER(LEN=256) :: var
    !WRITE(var,FMT="(A)") a
    PRINT*,"------DEBUG---"//TRIM(msg)//"=", a 
    RETURN
  end SUBROUTINE print_debug_msg_i
  ! ---------------------------------
  SUBROUTINE print_debug_msg_c( msg, a )
    ! write a dubug message to screen during running time
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: msg
    CHARACTER(LEN=*),INTENT(IN) :: a
    !CHARACTER(LEN=256) :: var
    !WRITE(var,FMT="(A)") a
    PRINT*,"------DEBUG---"//TRIM(msg)//"=", a 
    RETURN
  end SUBROUTINE print_debug_msg_c  


  ! ==========================================================================================

  SUBROUTINE report_error(content)
    implicit none
    CHARACTER(*),INTENT(IN)  :: content

    PRINT*,'Fatal error: ', content
    STOP
  end SUBROUTINE report_error

  ! ==========================================================================================

  SUBROUTINE percentage_display(n1, ntotal)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1, ntotal

    INTEGER :: k
    INTEGER,PARAMETER :: percents(9) = (/10, 20, 30, 40, 50, 60, 70, 80, 90/)
    INTEGER,PARAMETER :: percents2(4) = (/20, 40, 60, 80/)
    REAL(DPR),PARAMETER :: float_per(9) = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9/)
    REAL(DPR),PARAMETER :: float_per2(4) = (/0.2, 0.4, 0.6, 0.8/)

    IF( ntotal >= 100000) THEN
       DO k = 1, SIZE(float_per)
          IF( n1 == CEILING(ntotal * float_per(k)) ) THEN
             WRITE(*,'(I3, A2)') percents(k), ' %'
             RETURN
          END IF
       END DO
    ELSEIF( ntotal >= 50000 .AND. ntotal < 100000)  THEN
       DO k = 1, SIZE(float_per2)
          IF( n1 == CEILING(ntotal * float_per2(k)) ) THEN
             WRITE(*,'(I3, A2)') percents2(k), ' %'
             RETURN
          END IF
       END DO
    END IF

    IF( n1 == ntotal ) THEN
       WRITE(*,'(I3, A2)') 100, ' %'
    END IF

    RETURN
  end SUBROUTINE percentage_display

  ! ==========================================================================================

  !!===================================================================
  SUBROUTINE natural_log_spacing(ndata, minv, maxv, data)
    ! Note : maxv > minv > 0, ndata >1.
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: ndata
    REAL(DPR),INTENT(IN)  :: minv, maxv
    REAL(DPR),INTENT(INOUT)  :: data(ndata)
    REAL(DPR)  :: delta

    INTEGER  k

    delta = (LOG(maxv) - LOG(minv)) / (ndata -1)

    DO k = 1, ndata
       data(k) = EXP( LOG(minv) + delta * (k -1) )
    END DO

    RETURN
  END SUBROUTINE natural_log_spacing

  ! ---------------------------------

  SUBROUTINE linear_spacing(ndata, minv, maxv, data)
    ! Note : maxv > minv, ndata >1.
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: ndata
    REAL(DPR),INTENT(IN)  :: minv, maxv
    REAL(DPR),INTENT(INOUT)  :: data(ndata)
    REAL(DPR)  :: delta

    INTEGER  k

    delta = (maxv - minv) / (ndata -1)

    DO k = 1, ndata
       data(k) = minv + delta * (k -1)
    END DO

    RETURN

  END SUBROUTINE linear_spacing



  ! ---------------------------------

  SUBROUTINE meshing_grid_2D(n1, x1, x2, n2, y1, y2, array_x, array_y)
    ! Note : given two range pairs (e.g., [x0,x1] and [y0, y1]), and total amount of points in each direction of
    !  the two, generate a 2-D uniform mesh grid as a result of two 1-D arrays of X and Y, where:
    !  X: contains all x-values of the grid points,
    !  Y: contains all y-values of the grid points
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: n1, n2   ! >= 2
    REAL(DPR),INTENT(IN)  :: x1, x2, y1, y2
    REAL(DPR),INTENT(INOUT)  :: array_x( n1 * n2),  array_y( n1 * n2)
    REAL(DPR)  :: mx(n1), my(n2)

    INTEGER  k, k1, k2

    IF( (n1 <= 0) .OR. (n2 <= 0)) THEN
       PRINT*, ''
       PRINT*, '--In <meshing_grid_2D>, inputs of total amount of points are negative or zero ! '
       PRINT*, ' n1 == ', n1
       PRINT*, ' n2 == ', n2
       PRINT*, ' Running of program terminated by <STOP>'; STOP
    END IF


    CALL linear_spacing(n1, x1, x2, mx)
    CALL linear_spacing(n2, y1, y2, my)

    DO k = 1, n2
       k1 = 1 + (k-1) * n1
       k2 = k * n1

       array_x( k1 : k2 ) = mx    ! x- values
       array_y( k1 : k2 ) = my( k )  ! y- coordinate
    END DO

    RETURN

  END SUBROUTINE meshing_grid_2D


  ! ---------------------------------
  SUBROUTINE complex_dot_product( n, x1, x2, y )
    USE float_precision, ONLY: DPR
    IMPLICIT NONE
    INTEGER :: n, k
    COMPLEX(CDPR),INTENT(IN) :: x1(n), x2(n)
    COMPLEX(CDPR),INTENT(INOUT) :: y

    y = CMPLX(0.d0, 0.d0, KIND=DPR)

    DO k = 1, n
       y = y + x1(k) * x2(k)
    end DO
    RETURN
  end SUBROUTINE complex_dot_product




end MODULE FORTRAN_generic
