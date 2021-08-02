MODULE EM_data_postprocess
  ! Author: Jianbo Long, Jan, 2017
  USE float_precision, ONLY: DPR, CDPR
  USE error_cls
  USE FORTRAN_generic, ONLY: print_general_msg, print_debug_msg
  USE modelling_parameter, ONLY: modelling
  IMPLICIT NONE
  PRIVATE
  INTEGER,PARAMETER :: sub_len = 300
  INTEGER,PARAMETER :: wordlen = 256
  CHARACTER(LEN=sub_len),PROTECTED :: thisModule = "EM_data_postprocess"
  ! Public subroutines/functions
  PUBLIC :: data_processing_EM

  INTERFACE post_interpolation
     MODULE PROCEDURE post_interpolation_CM, post_interpolation_RE
  END INTERFACE post_interpolation
  
CONTAINS
  SUBROUTINE data_processing_EM ()
    USE file_rw, ONLY: write_numerical_result_to_files, write_to_files_Ron_MT2D
    USE EM_data_measurements, ONLY: nobs, xobs, yobs, zobs, get_synthetic_measurement_positions
    USE linear_system_equation_data, ONLY: EM_solution, gt_solution
    IMPLICIT NONE
    INTEGER :: num_rhs
    REAL(DPR),ALLOCATABLE :: fun_re(:)
    COMPLEX(CDPR),ALLOCATABLE :: fun_cm(:)
    ! double check the measurement locations ----
    IF( .NOT. ALLOCATED(xobs)) CALL get_synthetic_measurement_positions()

    IF( TRIM( ADJUSTL(modelling%DataType) ) == 'CM' ) THEN
       CALL print_general_msg("Writing MT2D sol for Ron's DD")
       CALL write_to_files_Ron_MT2D(EM_solution(:,1), EM_solution(:,2))
    end IF

    IF( nobs > 0 ) THEN
       IF( TRIM( ADJUSTL(modelling%DataType) ) == 'CM' ) THEN
          ALLOCATE(fun_cm(nobs) )
          NUM_RHS = SIZE(EM_solution, 2)
          CALL post_interpolation( nobs,xobs, zobs, fun_cm, EM_solution(:,1) )
          CALL write_numerical_result_to_files(nobs, xobs, yobs, zobs, fun_cm)
       ELSE
          ALLOCATE(fun_re(nobs) )
          NUM_RHS = SIZE(gt_solution, 2)
          CALL post_interpolation( nobs,xobs, zobs, fun_re, gt_solution(:,1) )
          CALL write_numerical_result_to_files(nobs, xobs, yobs, zobs, fun_re)
       END IF
    end IF
    RETURN
  end SUBROUTINE data_processing_EM

  SUBROUTINE post_interpolation_CM(nobs, xobs, zobs, func, solution)
    ! FE method over triangular meshes
    ! Inputs: # of data, 2D locations (x,z), and func component values at the dofs.
    ! Outputs: func at the measurements.
    USE FORTRAN_generic, ONLY: complex_dot_product
    USE mesh2d_discretization, ONLY: ele2node, ele2dof
    USE scalar_FE_kernels_2d, ONLY: which_triag_for_point,&
         get_basisWeights_linear_triag, get_shapeFunc_value_triag
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nobs
    REAL(DPR),INTENT(IN) :: xobs(nobs), zobs(nobs)
    COMPLEX(CDPR),INTENT(IN) :: solution(:)
    COMPLEX(CDPR),INTENT(INOUT):: func(nobs)
    REAL(DPR) :: BasisWeights(3,3)
    INTEGER :: i, k, nshape, whichele, ierr
    INTEGER :: localVerts(3)
    INTEGER,ALLOCATABLE :: localDofs(:)
    REAL(DPR) :: TempPoint(2), ele_area

    ! derivatives
    REAL(DPR),ALLOCATABLE :: Coe(:), CoeDx(:), CoeDz(:)
    COMPLEX(CDPR),ALLOCATABLE :: localCoe(:), localCoeDx(:), localCoeDz(:),&
         localFuncs(:)
    TYPE(error_type) :: error
    CHARACTER(LEN=wordlen) :: msg
    CHARACTER(LEN=sub_len) :: thisSubroutine = "post_interpolation_CM"

    CALL print_general_msg(TRIM(ADJUSTL(thisSubroutine)))
    !nfe = 1      ! order of FE
    IF(modelling%FEdegree == 1) THEN
       nshape = 3   ! local # of dofs
    ELSEIF(modelling%FEdegree == 2) THEN
       nshape = 6
    END IF
    ALLOCATE( localFuncs(nshape), localDofs(nshape) )
    ALLOCATE( Coe(nshape), CoeDx(nshape), CoeDz(nshape) )
    ALLOCATE( localCoe(nshape), localCoeDx(nshape), localCoeDz(nshape) )

    OBB: DO i = 1, nobs
       TempPoint(1) = xobs(i)
       TempPoint(2) = zobs(i)
       CALL which_triag_for_point( TempPoint, whichele, ierr )
       IF( ierr > 0 ) THEN  
          WRITE(msg,FMT='(A)') 'Possibly degenerated element'
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
          IF(error_check(error)) CALL error_report(error)
       END IF

       localVerts = ele2node( whichele, 1 : 3)
       localDofs = ele2dof( whichele, :)
       CALL get_basisWeights_linear_triag( localVerts, BasisWeights, ele_area)
       CALL get_shapeFunc_value_triag(nshape, BasisWeights, TempPoint, Coe, CoeDx, CoeDz)
       localCoe = CMPLX(Coe, 0.d0, KIND=8 )
       localCoeDx = CMPLX(CoeDx, 0.d0, KIND=8 )
       localCoeDz = CMPLX(CoeDz, 0.d0, KIND=8 )

       ! loop over the DOFs (or the basis functions)
       DO k = 1, nshape
          ! assemble local field values in vector form.
          localFuncs(k) = solution( localDofs(k) )  ! solution
       ENDDO

       CALL complex_dot_product( SIZE(localFuncs), localCoe, localFuncs, func(i) )

    end DO OBB
    RETURN
  end SUBROUTINE post_interpolation_CM

  SUBROUTINE post_interpolation_RE(nobs, xobs, zobs, func, solution)
    USE mesh2d_discretization, ONLY: ele2node, ele2dof
    USE scalar_FE_kernels_2d, ONLY: which_triag_for_point,&
         get_basisWeights_linear_triag, get_shapeFunc_value_triag
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nobs
    REAL(DPR),INTENT(IN) :: xobs(nobs), zobs(nobs)
    REAL(DPR),INTENT(IN) :: solution(:)
    REAL(DPR),INTENT(INOUT):: func(nobs)
    REAL(DPR) :: BasisWeights(3,3)
    INTEGER :: i, k, nshape, whichele, ierr
    INTEGER :: localVerts(3)
    INTEGER,ALLOCATABLE :: localDofs(:)
    REAL(DPR) :: TempPoint(2), ele_area

    ! derivatives
    REAL(DPR),ALLOCATABLE :: Coe(:), CoeDx(:), CoeDz(:), localFuncs(:)
    TYPE(error_type) :: error
    CHARACTER(LEN=wordlen) :: msg
    CHARACTER(LEN=sub_len) :: thisSubroutine = "post_interpolation_RE"

    CALL print_general_msg(TRIM(ADJUSTL(thisSubroutine)))
    IF(modelling%FEdegree == 1) THEN
       nshape = 3   ! local # of dofs
    ELSEIF(modelling%FEdegree == 2) THEN
       nshape = 6
    END IF
    ALLOCATE( localFuncs(nshape), localDofs(nshape) )
    ALLOCATE( Coe(nshape), CoeDx(nshape), CoeDz(nshape) )

    OBB: DO i = 1, nobs
       TempPoint(1) = xobs(i)
       TempPoint(2) = zobs(i)
       CALL which_triag_for_point( TempPoint, whichele, ierr )
       IF( ierr > 0 ) THEN  
          WRITE(msg,FMT='(A)') 'Possibly degenerated element'
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
          IF(error_check(error)) CALL error_report(error)
       END IF

       localVerts = ele2node( whichele, 1 : 3)
       localDofs = ele2dof( whichele, :)
       CALL get_basisWeights_linear_triag( localVerts, BasisWeights, ele_area)
       CALL get_shapeFunc_value_triag(nshape, BasisWeights, TempPoint, Coe, CoeDx, CoeDz)
       ! loop over the DOFs (or the basis functions)
       DO k = 1, nshape
          ! assemble local field values in vector form.
          localFuncs(k) = solution( localDofs(k) )  ! solution
       ENDDO

       func(i) = DOT_PRODUCT( Coe, localFuncs )
    end DO OBB
    RETURN
  end SUBROUTINE post_interpolation_RE

  ! -------------------------------------------------------
  ! -------------------------------------------------------


end MODULE EM_data_postprocess
