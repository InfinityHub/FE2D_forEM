MODULE num_algorithm_concatenate_differ_operator
  ! Jianbo Long, Mar, 2020
  USE modelling_parameter, ONLY: modelling
  IMPLICIT NONE
  PRIVATE
  ! public procedures
  PUBLIC :: differ_operator_algorithm_concatenate
CONTAINS
  SUBROUTINE differ_operator_algorithm_concatenate()
    ! to modify the data in the linear system module accessed by 'USE-CALL-procedure'
    USE linear_system_equation_data, ONLY: initialize_matrices_operator
    USE scalar_FE_kernels_2d, ONLY: get_submatrix_operator_linearFE_triag, &
         get_submatrix_operator_quadFE_triag
    IMPLICIT NONE
    CALL initialize_matrices_operator()

    IF(modelling%FEdegree == 1) THEN
       CALL get_submatrix_operator_linearFE_triag()
    ELSEIF(modelling%FEdegree == 2) THEN
       CALL get_submatrix_operator_quadFE_triag()
    END IF
    RETURN
  end SUBROUTINE differ_operator_algorithm_concatenate
end MODULE num_algorithm_concatenate_differ_operator
