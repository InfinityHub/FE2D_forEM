MODULE mesh_concatenate_module
  ! Author: Jianbo Long, Mar, 2020
  IMPLICIT NONE
  PRIVATE
  ! output
  INTEGER,PROTECTED :: ndof  
  PUBLIC :: mesh_concatenate
  PUBLIC :: ndof
CONTAINS
  ! -----------------------------------------------------------------------
  SUBROUTINE mesh_concatenate()
    USE mesh2d_discretization, ONLY: use_2D_triangular_mesh
    CALL use_2D_triangular_mesh(ndof)
    RETURN
  END SUBROUTINE Mesh_Concatenate
end MODULE Mesh_Concatenate_Module





