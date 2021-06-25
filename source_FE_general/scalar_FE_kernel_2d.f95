MODULE scalar_FE_kernels_2d
  ! Author: Jianbo Long, April, 2020
  USE float_precision, ONLY: DPR, CDPR
  USE error_cls
  USE strings_mod, ONLY: WORDLEN
  USE derived_data_module, ONLY: dofs

  USE FORTRAN_generic, ONLY: print_general_msg, print_debug_msg, add_to_list_integer
  USE linear_system_equation_data, ONLY: MatBandOne, MatBandLapPhi, &
       MatBand_div_resis_grad, MatBand_mu, Matband_dxdz

  IMPLICIT NONE
  PRIVATE
  INTEGER,PARAMETER :: sub_len = 300
  CHARACTER(LEN=sub_len),PROTECTED :: thisModule = "scalar_FE_kernels_2d"

  PUBLIC :: get_submatrix_operator_linearFE_triag, &
       get_submatrix_operator_quadFE_triag, &
       get_basisWeights_linear_triag, get_shapeFunc_value_triag, which_triag_for_point

CONTAINS
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  SUBROUTINE get_submatrix_operator_linearFE_triag()
    ! Get the non-zeros in each row of the final coefficient matrix in the
    ! nodal FE for each differential operators in the equations. Due to the nature of unstructured elements, the number
    ! of non-zeros in each row will not be a constant. An allocatable array
    ! is created for storing all non zeros.
    !
    ! Final outputs are in the banded matrix format
    !
    ! real-valued; triangular elements; linear-order
    !
    USE constants_module, ONLY: mu0
    USE mesh2d_discretization, ONLY: ele2node, p2node, p2ele, cellsigma, &
         set_cell_physical_property_2d
    USE mesh_concatenate_module, ONLY: ndof    ! # of dofs
    IMPLICIT NONE

    INTEGER  :: k, it, nelement, iele, nfe, nshape
    INTEGER,ALLOCATABLE :: Elemental_Dofs(:), mapIndex(:)
    INTEGER :: num_connected_dofs, num_connected_ele, testPointIndex
    INTEGER,ALLOCATABLE :: entireRowNodes(:)
    REAL(DPR) :: area_ele
    REAL(DPR),ALLOCATABLE :: BasisWeightsLinear(:,:), element_coef(:)
    !REAL(DPR),ALLOCATABLE :: vector_coef(:,:)
    REAL(DPR),ALLOCATABLE :: scattered_coef(:)

    !CHARACTER(LEN=300) :: subroutineTitle = "get_submatrix_operator_linearFE_triag"

    CALL print_general_msg('get submatrices of 2-D Diff Operators using nodal linear FE')

    nelement = SIZE( ele2node, 1 )
    CALL print_debug_msg(' # of total elements ', nelement)
    IF(.NOT. ALLOCATED(CellSigma) ) CALL set_cell_physical_property_2d('EM')

    nfe = 1;   nshape = 3
    ALLOCATE(Elemental_Dofs(nshape), mapIndex(nshape))
    ALLOCATE(BasisWeightsLinear( nshape, nshape ), element_coef(nshape))
    ! ALLOCATE(, vector_coef(nshape, 2))
    ! loop through all test functions (of interior dofs)
    PP:DO it = 1, ndof

       IF( dofs( it )%bud /= 1 )  THEN

          num_connected_dofs = SIZE( p2node(it)%nodes )
          num_connected_ele = SIZE( p2ele(it)%ele )

          IF( ALLOCATED(entireRowNodes) )  DEALLOCATE( entireRowNodes )
          ALLOCATE( entireRowNodes( num_connected_dofs)  )

          entireRowNodes = p2node(it)%nodes

          ALLOCATE( MatBandOne(it)%amat(num_connected_dofs) )
          ALLOCATE( MatBandOne(it)%jca(num_connected_dofs) )

          MatBandOne(it)%jca = entireRowNodes

          ALLOCATE( MatBandLapPhi(it)%amat(num_connected_dofs) )
          ALLOCATE( MatBand_div_resis_grad(it)%amat(num_connected_dofs) )
          ALLOCATE( MatBand_mu(it)%amat(num_connected_dofs) )

          IF( ALLOCATED(scattered_coef) )  DEALLOCATE( scattered_coef )
          ALLOCATE( scattered_coef( num_connected_dofs)  )

          ! loop through all connected elements for a node
          MatBandOne(it)%amat = 0.0
          MatBand_div_resis_grad(it)%amat = 0.0
          MatBand_mu(it)%amat = 0.0
          MatBandLapPhi(it)%amat = 0.0
          scattered_coef = 0.0

          Connected_cells: DO k = 1, num_connected_ele

             iele = p2ele(it)%ele( k )
             Elemental_Dofs = ele2node( iele, 1:3 )

             ! what is the index of the test point in the local dofs ?
             CALL get_mainNode_index( nshape, Elemental_Dofs, it, testPointIndex )
             ! mapping of the local dofs (3 for linear, 6 for quadratic) to the entire connnected dofs for a test point.
             CALL find_map_array1_to_array2( nshape, Elemental_Dofs, num_connected_dofs, entireRowNodes, mapIndex )
             ! get the weights for *linear* shape functions and element's area
             CALL get_basisWeights_linear_triag( Elemental_Dofs, BasisWeightsLinear, area_ele)

             ! (1) D = Laplacian
             CALL get_coef_Ni_dot_Laplacian_phi_triag( nfe, testPointIndex, area_ele, BasisWeightsLinear, nshape, element_coef )
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBandLapPhi(it)%amat = MatBandLapPhi(it)%amat + scattered_coef

             ! (2)D = Div(c* grad \phi) - \phi is a scalar function  (c is a physical property here): related to "D=Laplacian"
             element_coef = element_coef * (1.0/cellsigma( iele ))
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBand_div_resis_grad(it)%amat = MatBand_div_resis_grad(it)%amat + scattered_coef

             ! (3) D = 1 * mu * sigma
             CALL get_coef_Ni_dot_phi_triag( nfe, testPointIndex, area_ele, nshape, element_coef )
             ! physical properties (e.g., conductivity, permeability) are cell-based
             element_coef = element_coef * mu0 * cellsigma( iele )
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBandOne(it)%amat = MatBandOne(it)%amat + scattered_coef

             ! (4) D = 1 * mu
             CALL get_coef_Ni_dot_phi_triag( nfe, testPointIndex, area_ele, nshape, element_coef )
             ! physical properties (e.g., conductivity, permeability) are cell-based
             element_coef = element_coef * mu0
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBand_mu(it)%amat = MatBand_mu(it)%amat + scattered_coef

          end DO Connected_cells

       ENDIF
    end DO PP
    RETURN
  end SUBROUTINE get_submatrix_operator_linearFE_triag
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  SUBROUTINE get_submatrix_operator_quadFE_triag()
    ! 
    ! real-valued; triangular elements; quadratic-order
    !
    USE constants_module, ONLY: mu0
    USE modelling_parameter, ONLY: modelling
    USE mesh2d_discretization, ONLY: ele2node, p2node, p2ele, edge2node, cellsigma, &
         ele2dof, eg2eg, eg2ele, set_cell_physical_property_2d
    USE mesh_concatenate_module, ONLY: ndof    ! # of dofs
    IMPLICIT NONE

    INTEGER  :: k, it, nelement, iele, nfe, nshape, nedge, nvert, id_loc
    INTEGER,ALLOCATABLE :: list_dofs(:), list_ele(:)
    INTEGER :: num_connected_dofs, num_connected_ele, testPointIndex, Elemental_verts(3)
    INTEGER,ALLOCATABLE :: Elemental_dofs(:), mapIndex(:)
    REAL(DPR) :: area_ele, BasisWeightsLinear(3,3)
    REAL(DPR),ALLOCATABLE :: element_coef(:)
    REAL(DPR),ALLOCATABLE :: scattered_coef(:)
    INTEGER,ALLOCATABLE :: extra_edges(:), extra_verts(:)

    CALL print_general_msg('get submatrices of 2-D Diff Operators using nodal quad FE')

    nelement = SIZE( ele2node, 1 )
    nedge = SIZE(edge2node, 1)
    nvert = ndof - nedge    ! NO. of vertices
    CALL print_debug_msg(' # of total elements ', nelement)
    IF(modelling%DataType == "CM") THEN
       IF(.NOT. ALLOCATED(CellSigma) ) CALL set_cell_physical_property_2d('EM')
    END IF

    nfe = 0; nshape = 0
    IF(modelling%FEdegree == 1) THEN
       nfe = 1;   nshape = 3
    ELSEIF(modelling%FEdegree == 2) THEN
       nfe = 2;   nshape = 6
    end IF

    ALLOCATE(Elemental_dofs(nshape), mapIndex(nshape))
    ALLOCATE(element_coef(nshape))

    ! loop through all test functions
    PP:DO it = 1, ndof
       ! boundary dofs, just cycle the loop
       IF( dofs(it)%bud == 1) CYCLE PP

       IF( it <= nvert) THEN
          IF( ALLOCATED(extra_edges) )  DEALLOCATE( extra_edges )
          CALL connected_edges_for_vert_search( it, extra_edges )
          num_connected_dofs = SIZE( p2node(it)%nodes ) + SIZE(extra_edges)
          num_connected_ele = SIZE( p2ele(it)%ele )
          IF( ALLOCATED(list_dofs) )  DEALLOCATE( list_dofs )
          ALLOCATE(list_dofs( num_connected_dofs ))
          list_dofs = [p2node(it)%nodes, extra_edges + nvert] ! global numbering
          IF( ALLOCATED(list_ele) )  DEALLOCATE( list_ele )
          ALLOCATE(list_ele(num_connected_ele))
          list_ele = p2ele(it)%ele
       ELSE
          id_loc = it - nvert
          IF( ALLOCATED(extra_verts) )  DEALLOCATE( extra_verts )
          CALL connected_verts_for_edge_search( id_loc, extra_verts )
          num_connected_dofs = SIZE( eg2eg(id_loc)%edge ) + SIZE(extra_verts)
          num_connected_ele = SIZE( eg2ele(id_loc)%ele )
          IF( ALLOCATED(list_dofs) )  DEALLOCATE( list_dofs )
          ALLOCATE(list_dofs(num_connected_dofs))
          list_dofs = [eg2eg(id_loc)%edge + nvert,  extra_verts] ! global numbering
          IF( ALLOCATED(list_ele) )  DEALLOCATE( list_ele )
          ALLOCATE(list_ele(num_connected_ele))
          list_ele = eg2ele(id_loc)%ele
       END IF

       ALLOCATE( MatBandOne(it)%amat(num_connected_dofs) )
       ALLOCATE( MatBandOne(it)%jca(num_connected_dofs) )

       MatBandOne(it)%jca = list_dofs

       IF(modelling%DataType == "CM") THEN
          ALLOCATE( MatBandLapPhi(it)%amat(num_connected_dofs) )
          ALLOCATE( MatBand_div_resis_grad(it)%amat(num_connected_dofs) )
          ALLOCATE( MatBand_mu(it)%amat(num_connected_dofs) )
       ELSE
          ALLOCATE( MatBand_dxdz(it)%amat(num_connected_dofs) )
       END IF

       IF( ALLOCATED(scattered_coef) )  DEALLOCATE( scattered_coef )
       ALLOCATE( scattered_coef( num_connected_dofs)  )

       ! loop through all connected elements for a node
       IF(modelling%DataType == "CM") THEN
          MatBandOne(it)%amat = 0.0
          MatBand_div_resis_grad(it)%amat = 0.0
          MatBand_mu(it)%amat = 0.0
          MatBandLapPhi(it)%amat = 0.0
       ELSE
          MatBand_dxdz(it)%amat = 0.0
       END IF
       scattered_coef = 0.0

       Connected_cells: DO k = 1, num_connected_ele

          iele = list_ele( k )
          ! elemental vertices are required for computing linear shape functions
          Elemental_verts = ele2node( iele, 1:3 )
          ! now global locations of elemental dofs
          Elemental_dofs = ele2dof( iele, :)

          ! what is the index of the test point in the local nodes (dofs) ?
          CALL get_mainNode_index( nshape, Elemental_dofs, it, testPointIndex )
          ! mapping of the local dof (3 for linear, 6 for quadratic) to the entire connnected dofs for a test point.
          CALL find_map_array1_to_array2( nshape, Elemental_dofs, num_connected_dofs, list_dofs, mapIndex )
          ! get the weights for *linear* shape functions (regardless of high-order FE) and element's area
          CALL get_basisWeights_linear_triag( Elemental_verts, BasisWeightsLinear, area_ele)

          IF(modelling%DataType == "CM") THEN
             ! (1) D = Laplacian
             CALL get_coef_Ni_dot_Laplacian_phi_triag( nfe, testPointIndex, area_ele, BasisWeightsLinear, nshape, element_coef )
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBandLapPhi(it)%amat = MatBandLapPhi(it)%amat + scattered_coef

             ! (2)D = Div(c* grad \phi) - \phi is a scalar function  (c is a physical property here)
             element_coef = element_coef * (1.0/cellsigma( iele ))
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBand_div_resis_grad(it)%amat = MatBand_div_resis_grad(it)%amat + scattered_coef

             ! (3) D = 1 * mu * sigma
             CALL get_coef_Ni_dot_phi_triag( nfe, testPointIndex, area_ele, nshape, element_coef )
             ! physical properties (e.g., conductivity, permeability) are cell-based
             element_coef = element_coef * mu0 * cellsigma( iele )
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBandOne(it)%amat = MatBandOne(it)%amat + scattered_coef

             ! (4) D = 1 * mu
             CALL get_coef_Ni_dot_phi_triag( nfe, testPointIndex, area_ele, nshape, element_coef )
             ! physical properties (e.g., conductivity, permeability) are cell-based
             element_coef = element_coef * mu0
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBand_mu(it)%amat = MatBand_mu(it)%amat + scattered_coef

          ELSE
             ! (5) D = dxdz (2nd order);  no physical property
             CALL get_coef_Ni_dot_dxdz_phi_triag( nfe, testPointIndex, area_ele, BasisWeightsLinear, nshape, element_coef )
!!$             IF( COUNT(ABS(element_coef) > 1.E-10 ) <= 0 ) THEN
!!$                PRINT*, "too few non-zero numerics !", it, k, num_connected_ele
!!$                PRINT*, "testPointIndex ", testPointIndex
!!$                !PRINT*, "amat is: ", element_coef
!!$                PRINT*, "Linear Basis L1:", BasisWeightsLinear(1,:)
!!$                PRINT*, "Linear Basis L2:", BasisWeightsLinear(2,:)
!!$                PRINT*, "Linear Basis L3:", BasisWeightsLinear(3,:)
!!$                !READ*, id_loc
!!$             end IF
             CALL mapping_array1_to_array2( nshape, element_coef, mapIndex, num_connected_dofs, scattered_coef )
             MatBand_dxdz(it)%amat = MatBand_dxdz(it)%amat + scattered_coef
          END IF
       end DO Connected_cells

    end DO PP
    RETURN
  end SUBROUTINE get_submatrix_operator_quadFE_triag




  !---------------------------------------------------------------
  !---------------------------------------------------------------  

  SUBROUTINE get_coef_Ni_dot_Laplacian_phi_triag( nfe, testPointIndex, area, linearCoef, nshape, int_coeff )
    ! Get coefficients in the Integral: Ni * nabla^2 Phi dV over a triangle.
    !  the scalar phi can be any scalar function (e.g., Ax,Ay,Az & \phi in the A-\phi scheme)
    !
    ! \int Ni * nabla^2 Phi dV = - \int (grad Ni) * (grad phi) dV  when the surface integration
    !  is cancelled out.
    !
    ! This integral is the special case of \int Ni * div(c* grad phi) dV where c is a scalar function
    !    (the spcial case c=1.)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: testPointIndex    ! local index of the main node associated with Ni
    INTEGER,INTENT(IN) :: nfe
    INTEGER,INTENT(IN) :: nshape       ! number of shape functions within a cell
    REAL(DPR),INTENT(IN) :: area, linearCoef(3,3)
    REAL(DPR),INTENT(OUT) :: int_coeff(nshape) ! local coefficients resulted from the volume integration
    INTEGER :: k

    DO k = 1, nshape
       CALL get_SI_gradNi_gradNk_triag( nfe, testPointIndex, k, area, linearCoef, int_coeff(k) )
    END DO

    ! the minus sign goes here to make the result be consistent with the subroutine's name
    int_coeff = int_coeff * (-1)
    RETURN
  end SUBROUTINE get_coef_Ni_dot_Laplacian_phi_triag


  !---------------------------------------------------------------
  !---------------------------------------------------------------  

  SUBROUTINE get_coef_Ni_dot_grad_phi_triag( nfe, testPointIndex, area, linearCoef, nshape, int_coeff )
    ! Get coefficients in the Integral: Ni * \nabla Phi dV over a triangle.
    !  the scalar phi can be any scalar function (e.g., Ax,Ay,Az & \phi in the A-\phi scheme)
    !
    ! \int Ni * \nabla Phi dV = \SUM(\int Ni * \nabla Nk) dV 
    !
    ! This subroutine is related to \int A \cdot grad phi dV for any vector A and scalar phi.
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: testPointIndex   ! local index of the test point associated with Ni
    INTEGER,INTENT(IN) :: nfe            ! degree of FE basis functions
    INTEGER,INTENT(IN) :: nshape          ! number of shape functions within a cell
    REAL(DPR),INTENT(IN) :: area, linearCoef(3,3)
    REAL(DPR),INTENT(OUT) :: int_coeff(nshape, 2) ! local coefficients resulting from the area integration

    INTEGER :: k
    REAL(DPR) :: temp_coef(2)

    DO k = 1, nshape
       CALL get_SI_Ni_gradNk_triag( nfe, testPointIndex, k, area, linearCoef, temp_coef )
       int_coeff(k, :) = temp_coef
    END DO
    RETURN
  end SUBROUTINE get_coef_Ni_dot_grad_phi_triag
  !---------------------------------------------------------------
  !---------------------------------------------------------------  

  SUBROUTINE get_coef_Ni_dot_DivA_triag( nfe, testPointIndex, area, linearCoef, nshape, int_coeff )
    ! Get coefficients in the Integral: Ni * (\div A) dA over a triangle.
    !  A is a vector.
    !
    ! \int Ni * (\div A) dV =  "Surface_integration" - \int A \cdot \nabla Ni dV
    !   = - {\int Ax * (grad Ni)_x dV,  \int Az * (grad Ni)_z dV}
    !
    !  where (grad Ni)_x = dNi/dx
    ! (Surface_integration is assumed to be zero)
    !
    ! And for each component of the above:
    !     \int f * grad Ni dV = \sum^{ns}_{k=1} {\int Nk * grad Ni dV}
    ! by expanding f in terms of <ns> shape functions.
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: testPointIndex    ! local index of the test point associated with Ni
    INTEGER,INTENT(IN) :: nfe             ! degree of FE basis functions
    INTEGER,INTENT(IN) :: nshape           ! number of shape functions within a cell
    REAL(DPR),INTENT(IN) :: area, linearCoef(3,3)
    REAL(DPR),INTENT(OUT) :: int_coeff(nshape, 2) ! local coefficients resulting from the 2D integration

    INTEGER :: k
    REAL(DPR) :: temp_coef(2)

    DO k = 1, nshape
       ! note below, i and k are switched, so that it can return Nk_gradNi, rather Ni_gradNk
       CALL get_SI_Ni_gradNk_triag( nfe, k, testPointIndex, area, linearCoef, temp_coef )
       int_coeff(k, 1:2) = temp_coef
    END DO

    int_coeff = int_coeff * (-1)

    RETURN
  end SUBROUTINE get_coef_Ni_dot_DivA_triag


  !---------------------------------------------------------------
  !---------------------------------------------------------------

  SUBROUTINE get_coef_Ni_dot_phi_triag( nfe, testPointIndex, area, nshape, int_coeff )
    ! Get coefficients in the Integral: Ni * Phi dA over a triangle.
    !  the scalar phi can be any scalar function.
    !
    !  int Ni * Phi dV = [c1 c2 c3] * [phi1 phi2 phi3]^T
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nfe
    INTEGER,INTENT(IN) :: testPointIndex    ! local index for the test point (<= nshape)
    INTEGER,INTENT(IN) :: nshape           ! number of shape functions within a cell
    REAL(DPR),INTENT(IN) :: area
    REAL(DPR),INTENT(OUT) :: int_coeff(nshape) ! local coefficients resulting from integration
    INTEGER :: k

    ! nshape ==
    !  3  for linear-order (if nfe =1)
    !  6 for quadratic-order (if nfe =2)
    DO k = 1, nshape
       CALL get_SI_Ni_Nk_triag( nfe, testPointIndex, k, area, int_coeff(k) )
    END DO
    RETURN
  end SUBROUTINE get_coef_Ni_dot_phi_triag
  !---------------------------------------------------------------
  !---------------------------------------------------------------

  SUBROUTINE get_coef_Ni_dot_dxdz_phi_triag( nfe, testPointIndex, area, linearCoef, nshape, int_coeff )
    ! Get coefficients in the Integral: Ni * d^2(Phi)/dxdz dA over a triangle.
    !  the scalar phi can be any scalar function.
    !
    !  int Ni * Phi dV = [c1 c2 c3] * [phi1 phi2 phi3]^T, for example.
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nfe
    INTEGER,INTENT(IN) :: testPointIndex    ! local index for the test point (<= nshape)
    INTEGER,INTENT(IN) :: nshape           ! number of shape functions within a cell
    REAL(DPR),INTENT(IN) :: area, linearCoef(3,3)
    REAL(DPR),INTENT(OUT) :: int_coeff(nshape) ! local coefficients resulting from integration
    INTEGER :: k

    ! nshape ==
    !  3  for linear-order (if nfe =1)
    !  6 for quadratic-order (if nfe =2)
    DO k = 1, nshape
       CALL get_SI_Ni_dxdz_Nk_triag( nfe, testPointIndex, k, area, linearCoef, int_coeff(k) )
    END DO
    RETURN
  end SUBROUTINE get_coef_Ni_dot_dxdz_phi_triag

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  SUBROUTINE get_mainNode_index( np, AllNodes, mainNode, LocalIndex )
    ! Get the position of the main node (associated with the test function Ni) in
    !  the local array of vertices.
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: np
    INTEGER,INTENT(IN) :: AllNodes(np)      ! global nodal # of points within a considered element
    INTEGER,INTENT(IN) :: mainNode          ! global nodal # associated with Ni
    INTEGER,INTENT(OUT) :: LocalIndex   !  local index of the main node
    INTEGER :: k, whichIndex
    LOGICAL :: foundnode

    foundnode = .FALSE.

    FindNode: DO k = 1, np
       IF( mainNode == AllNodes(k) )  THEN
          ! The local index of the test node in the array. E.g., 3(the position) for 54 in /21, 3, 54, 78/
          whichIndex = k
          foundnode = .TRUE.
          EXIT FindNode
       END IF
    END DO FindNode

    IF( foundnode .EQV. .FALSE.) THEN
       PRINT*,''
       PRINT*,'Error in finding the local index for the main node ... Terminated !'
       PRINT*,'The main node inputed is: ', mainNode
       PRINT*, 'All local nodes inputed are: ',  AllNodes
       PRINT*, 'Program terminated by <STOP>'
       STOP
    END IF

    LocalIndex = whichIndex
    RETURN
  end SUBROUTINE get_mainNode_index

  !---------------------------------------------------------------
  !---------------------------------------------------------------
  SUBROUTINE get_SI_Ni_triag( nfe, i, area, integral_out )
    ! Integral: Ni dV , over a triangle, with
    !     Ni as scalar shape function for the test point/node
    !
    ! nfe: the polynomial degree (1: linear; 2: quadratic, 3: cubic etc)
    ! Based on the formulae for high-order FE shape functions in Jin's book, 3rd edition
    !
    ! Assuming scalar basis functions in each triangle
    !
    ! Jianbo Long, April, 2020
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nfe  ! degree of the polynomials (basis functions)
    INTEGER,INTENT(IN) :: i   ! local index of shape function within a triangle
    REAL(DPR),INTENT(IN) :: area
    REAL(DPR),INTENT(OUT) :: integral_out

    ! intermediate arguments
    REAL(DPR) :: t1, t2

    SELECT CASE( nfe )
    CASE(1)
       CALL get_linear_SI_N123_triag( 1, 0, 0, area, integral_out )

    CASE(2)
       ! The shape function is a Quadratic order shape function
       ! There are 2 situations regarding i
       IF( i >= 1 .AND. i <= 3 ) THEN
          ! situation 1: vertex nodes
          CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t1 )
          CALL get_incomplete_linear_SI_N123_triag( 1, 0, 0, t2 )
          integral_out = (2 * t1 - t2) * (2.0 * area)
       ELSEIF( i >= 4 .AND. i <= 6 ) THEN
          CALL get_linear_SI_N123_triag( 1, 1, 0, area, t1 )
          integral_out = t1 * 4.0
       END IF
    end SELECT
    RETURN
  end SUBROUTINE get_SI_Ni_triag
  !---------------------------------------------------------------
  SUBROUTINE get_SI_Ni_dxdz_Nk_triag( nfe, i, k, area, linearCoef, integral_out )
    ! Integral: Ni * d^2(Nk)/dxdz dV , over a triangle, with
    !     Ni as scalar shape function for the test point/node;
    !     d^2(Nk)/dxdz is the 2nd-order derivative of shape func Nk, and
    !     will be a non-zero constant if nfe=2; and 0 if nfe=1. 
    !
    ! nfe: the polynomial degree (1: linear; 2: quadratic, 3: cubic etc)
    ! This is used for unusual PDEs
    !
    ! Assuming scalar basis functions in each triangle
    ! Jianbo Long, August, 2020
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nfe  ! degree of the polynomials (basis functions)
    INTEGER,INTENT(IN) :: i, k   ! local indexes of shape functions (Ni, Nk) within a triangle
    REAL(DPR),INTENT(IN) :: area, linearCoef(3,3)
    REAL(DPR),INTENT(OUT) :: integral_out

    ! intermediate arguments
    INTEGER :: p, q
    REAL(DPR) :: t1, t2, dev

    SELECT CASE( nfe )
    CASE(1)
       integral_out = 0.0
    CASE(2)
       ! The shape function is a Quadratic order shape function
       ! Result of 2nd-order derivative of a quadratic shape func is a constant
       IF( k >= 1 .AND. k <= 3 ) THEN
          ! Type-1 shape func Nk: vertex-type
          dev = 4.0 * linearCoef(k,2) * linearCoef(k,3)
       ELSEIF( k >= 4 .AND. k <= 6 ) THEN
          ! Type-2 shape func Nk: edge-type
          CALL get_two_vertices_from_mid_edge_point(k, p, q)  ! Nk(x,z) = 4*Lp *Lq
          dev = 4.0 * ( linearCoef(p,2) * linearCoef(q,3) + linearCoef(p,3) * linearCoef(q,2) )
       END IF

       ! There are also 2 situations regarding Ni
       IF( i >= 1 .AND. i <= 3 ) THEN
          ! vertex-type Ni and Nk
          CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t1 )
          CALL get_incomplete_linear_SI_N123_triag( 1, 0, 0, t2 )

          integral_out = (2 * t1 - t2) * (2.0 * area) * dev

          IF(ABS(integral_out) < 1.E-10) THEN
             PRINT*, "t1:", t1
             PRINT*, "t2:", t2
             PRINT*, "(2 * t1 - t2) is: ", (2 * t1 - t2)
             PRINT*, "dev and (i,k) is: ", dev, i,k
             !READ*, id_loc
          end IF
       ELSEIF( i >= 4 .AND. i <= 6 ) THEN
          CALL get_linear_SI_N123_triag( 1, 1, 0, area, t1 )
          integral_out = t1 * 4.0 * dev
       END IF
    end SELECT
    RETURN
  end SUBROUTINE get_SI_Ni_dxdz_Nk_triag

  ! -----------------------------------------------------------

  SUBROUTINE get_SI_Ni_Nk_triag( nfe, i, k, area, integral_out )
    ! Integral: Ni * Nk dV , over a triangle, with
    !     Ni as scalar shape function for the test point/node,
    !     Nk as the shape function for the second point/node (maybe the same as the test point if k==i)
    !
    ! nfe: the polynomial degree (1: linear; 2: quadratic, 3: cubic etc)
    ! Based on the formulae for high-order FE shape functions in Jin's book, 3rd edition
    !
    ! Assuming scalar basis functions in each triangle
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nfe  ! degree of the polynomials (basis functions)
    INTEGER,INTENT(IN) :: i, k  ! local indices of shape functions within a triangle
    REAL(DPR),INTENT(IN) :: area
    REAL(DPR),INTENT(OUT) :: integral_out

    ! intermediate arguments
    REAL(DPR) :: t1, t2, t3, t4
    INTEGER :: m, n, p, q
    TYPE(error_type) :: error
    CHARACTER(LEN=wordlen) :: msg
    CHARACTER(LEN=sub_len) :: thisSubroutine = "get_SI_Ni_Nk_triag"
    SELECT CASE( nfe )
    CASE(1)
       ! the shape functions Ni,Nk are Linear shape functions
       IF( i >= 1 .AND. i <= 3 .AND. k >= 1 .AND. k <= 3 ) THEN
          IF( i /= k ) THEN
             CALL get_linear_SI_N123_triag( 1, 1, 0, area, integral_out )
          ELSE
             CALL get_linear_SI_N123_triag( 2, 0, 0, area, integral_out )
          ENDIF
       ELSE
          WRITE(msg,FMT='(A,A,A)') ' Invalid input of i and k (both should be >=1 and <=3) when <nfe = 1>: ', i,k
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
          IF(error_check(error)) CALL error_report(error)
       END IF

    CASE(2)
       ! The shape functions are Quadratic order shape functions
       ! There are 4 situations regarding i & k
       IF( i >= 1 .AND. i <= 3 .AND. k >= 1 .AND. k <= 3 ) THEN
          ! situation 1: both points are vertex nodes
          IF( i /= k ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 2, 2, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t2 )
             t3 = t2
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t4 )
          ELSE
             CALL get_incomplete_linear_SI_N123_triag( 4, 0, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 3, 0, 0, t2 )
             t3 = t2
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t4 )
          ENDIF

          integral_out = (4.0 * t1 - 2 * t2 - 2 * t3 + t4 ) * (2.0 * area)

       ELSEIF( i >= 1 .AND. i <= 3 .AND. k >= 4 .AND. k <= 6 ) THEN
          ! situation 2: test node is a vertex, the other is a mid-edge node
          CALL get_two_vertices_from_mid_edge_point(k, m, n)
          ! m /= n, but one of [m,n] can be equal to i, or none of them is equal to i
          IF( (m == i) .OR. (n == i) ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 3, 1, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t2 )
          ELSE
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 1, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 1, t2 )
          ENDIF

          integral_out = 4.0 * (2 * t1 - t2 ) * (2.0 * area)

       ELSEIF( i >= 4 .AND. i <= 6 .AND. k >= 1 .AND. k <= 3 ) THEN
          ! situation 3: test node is a mid-edge position, the other is a vertex
          CALL get_two_vertices_from_mid_edge_point(i, m, n)
          IF( (m == k) .OR. (n == k) ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 3, 1, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t2 )
          ELSE
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 1, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 1, t2 )
          ENDIF

          integral_out = 4.0 * (2 * t1 - t2 ) * (2.0 * area)

       ELSEIF( i >= 4 .AND. i <= 6 .AND. k >= 4 .AND. k <= 6 ) THEN
          ! situation 4: both nodes are mid-edge node
          CALL get_two_vertices_from_mid_edge_point(i, m, n)  ! Ni(x,y,z)
          CALL get_two_vertices_from_mid_edge_point(k, p, q)  ! Nk(x,y,z)
          ! The two edges can be : (1) the same; or (2) connected
          IF( i == k ) THEN
             ! two edges are the same: m == p, and n == q
             CALL get_incomplete_linear_SI_N123_triag( 2, 2, 0, t1 )  ! Lm*Ln*Lp*Lq
          ELSE
             ! two different edges (are connected)
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 1, t1 )
          ENDIF

          integral_out = 16.0 * t1 * (2.0 * area)

       ELSE
          WRITE(msg,FMT='(A,A,A)') ' Invalid input of i and k (both should be >=1 and <=6) when <nfe = 2>: ', i,k
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
          IF(error_check(error)) CALL error_report(error)
       END IF

    CASE DEFAULT
       WRITE(msg,FMT='(A,A)') ' Invalid input of <nfe> (should be 1 <= n <= 2 by far: ', nfe
       CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
       IF(error_check(error)) CALL error_report(error)
    end SELECT

    RETURN   
  end SUBROUTINE get_SI_Ni_Nk_triag

  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------

  SUBROUTINE get_SI_gradNi_gradNk_triag( nfe, i, k, area, linearCoef, integral_out )
    ! Integral: grad(Ni) * grad(Nk) dV , over a triangle, with
    !     Ni as scalar shape function for the test point/node,
    !     Nk as the shape function for the second point/node (maybe the same as the test point if k==i)
    !     grad(Ni) : the gradient of Ni
    !     grad(Nk) : the gradient of Nk
    !
    ! nfe: the polynomial degree (1: linear; 2: quadratic, 3: cubic etc)
    ! (Based on the formulae for high-order FE shape functions in Jin's book, 3rd edition)
    !
    ! linearCoef: the coefficients of 3 linear shape functions, pre-calculated
    !
    ! Jianbo Long, April, 2020
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nfe  ! degree of the polynomials (basis functions)
    INTEGER,INTENT(IN) :: i, k  ! local indices of shape functions within a triangle
    REAL(DPR),INTENT(IN) :: area, linearCoef(3,3)
    REAL(DPR),INTENT(OUT) :: integral_out

    ! intermediate arguments
    REAL(DPR) :: t1, t2, t3, t4, grad1, grad2, grad3, grad4
    INTEGER :: m, n, p, q
    TYPE(error_type) :: error
    CHARACTER(LEN=wordlen) :: msg
    CHARACTER(LEN=sub_len) :: thisSubroutine = "get_SI_gradNi_gradNk_triag"

    SELECT CASE( nfe )
    CASE(1)
       ! the shape functions Ni,Nk are Linear shape functions
       IF( (i >= 1) .AND. (i <= 3) .AND. (k >= 1) .AND. (k <= 3) ) THEN
          ! the gradient of Ni is contant
          integral_out = DOT_PRODUCT( linearCoef(i, 2:3), linearCoef(k, 2:3) ) * area

       ELSE
          WRITE(msg,FMT='(A,A,A)') ' Invalid input of i and k (both should be >=1 and <=3) when <nfe = 1>: ',i,k
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
          IF(error_check(error)) CALL error_report(error)
       END IF

    CASE(2)       
       ! The shape functions are Quadratic order shape functions
       ! There are 4 situations regarding i & k
       IF( (i >= 1) .AND. (i <= 3) .AND. (k >= 1) .AND. (k <= 3) ) THEN
          ! situation 1: both points are vertex nodes
          grad1 = DOT_PRODUCT( linearCoef(i, 2:3), linearCoef(k, 2:3) )

          IF( i /= k ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t1 )
          ELSE
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t1 )
          ENDIF

          CALL get_incomplete_linear_SI_N123_triag( 1, 0, 0, t2 )
          t3 = t2

          t4 = area

          integral_out = (16 * t1 - 4 * t2 - 4 * t3 ) * (2.0 * area) + t4
          integral_out = integral_out * grad1

       ELSEIF( (i >= 1) .AND. (i <= 3) .AND. (k >= 4) .AND. (k <= 6) ) THEN
          ! situation 2: test node is a vertex, the other is a mid-edge node
          CALL get_two_vertices_from_mid_edge_point(k, m, n)

          grad1 = DOT_PRODUCT( linearCoef(i, 2:3), linearCoef(m, 2:3) )
          grad2 = DOT_PRODUCT( linearCoef(i, 2:3), linearCoef(n, 2:3) )

          ! m /= n, but one of (m,n) can be equal to i, or none of them is equal to i
          IF( (m == i) ) THEN
             ! for Li * Ln (t1) and Li * Lm (t3)
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t3 )
          ELSEIF( n == i ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t3 )
          ELSEIF( (m/=i) .AND. (n/=i) )  THEN
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t3 )
          ENDIF

          CALL get_incomplete_linear_SI_N123_triag( 1, 0, 0, t2 ) ! Ln
          t4 = t2       ! Lm

          integral_out = 4.0 * ( (4 * t1 - t2) * grad1 + (4 *t3 - t4) * grad2 ) * (2.0 * area)

       ELSEIF( (i >= 4) .AND. (i <= 6) .AND. (k >= 1) .AND. (k <= 3) ) THEN
          ! situation 3: test node is a mid-edge position, the other is a vertex
          CALL get_two_vertices_from_mid_edge_point(i, m, n)

          grad1 = DOT_PRODUCT( linearCoef(m, 2:3), linearCoef(k, 2:3) )
          grad2 = DOT_PRODUCT( linearCoef(n, 2:3), linearCoef(k, 2:3) )

          IF( (m == k) ) THEN
             ! for Ln * Lk (t1) and Lm * Lk (t3)
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t3 )
          ELSEIF( n == k ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t3 )
          ELSEIF( (m/=k) .AND. (n/=k) ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t3 )             
          ENDIF

          CALL get_incomplete_linear_SI_N123_triag( 1, 0, 0, t2 )  ! Ln
          t4 = t2

          integral_out = 4.0 * ( (4 * t1 - t2) * grad1 + (4 *t3 - t4) * grad2 ) * (2.0 * area)

       ELSEIF( (i >= 4) .AND. (i <= 6) .AND. (k >= 4) .AND. (k <= 6) ) THEN
          ! situation 4: both nodes are mid-edge node
          CALL get_two_vertices_from_mid_edge_point(i, m, n)  ! Ni(x,y,z), m/=n
          CALL get_two_vertices_from_mid_edge_point(k, p, q)  ! Nk(x,y,z), p/=q
          ! The two edges can be : (1) the same; or (2) connected;

          grad1 = DOT_PRODUCT( linearCoef(m, 2:3), linearCoef(p, 2:3) )
          grad2 = DOT_PRODUCT( linearCoef(m, 2:3), linearCoef(q, 2:3) )
          grad3 = DOT_PRODUCT( linearCoef(n, 2:3), linearCoef(p, 2:3) )
          grad4 = DOT_PRODUCT( linearCoef(n, 2:3), linearCoef(q, 2:3) )

          nes:IF( i == k ) THEN
             ! two edges are the same: m == p, and n == q
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t1 ) ! Ln*Lq
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t2 ) ! Ln*Lp
             t3 = t2  ! Lm * Lq
             t4 = t1  ! Lm * Lp
          ELSE
             ! two edges are connected, creating 4 scenarios further, remember: m/=n, p/=q
             IF( m == p ) THEN
                ! we have: m(p), n, q three different vertices
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t1 ) ! Ln*Lq
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t2 ) ! Ln*Lp
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t3 ) ! Lm*Lq
                CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t4 ) ! Lm*Lp
             ELSEIF( m == q ) THEN
                ! NOW: m(q), n, p three different vertices
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t1 ) ! Ln*Lq
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t2 ) ! Ln*Lp
                CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t3 ) ! Lm*Lq
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t4 ) ! Lm*Lp
             ELSEIF( (m/=p) .AND. (m/=q) .AND. (n==p) ) THEN
                ! NOW: m, p, q three different vertices, n is one of p,q
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t1 ) ! Ln*Lq
                CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t2 ) ! Ln*Lp
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t3 ) ! Lm*Lq
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t4 ) ! Lm*Lp
             ELSEIF( n == q ) THEN
                ! NOW: m/=p,  m/=q, n/=p
                CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t1 ) ! Ln*Lq
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t2 ) ! Ln*Lp
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t3 ) ! Lm*Lq
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t4 ) ! Lm*Lp                   
             ENDIF
          END IF nes

          integral_out = 16.0 * (t1 * grad1 + t2 * grad2 + t3 * grad3 + t4 * grad4) * (2.0 * area)
       ELSE
          WRITE(msg,FMT='(A,A,A)') ' Invalid input of i and k (both should be >=1 and <=6) when <nfe = 2>: ',i,k
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
          IF(error_check(error)) CALL error_report(error)
       END IF
    CASE DEFAULT
       WRITE(msg,FMT='(A,A)') ' Invalid input of <nfe> (should be >=1 and <=2) by far: ', nfe
       CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
       IF(error_check(error)) CALL error_report(error)
    end SELECT

    RETURN
  end SUBROUTINE get_SI_gradNi_gradNk_triag

  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------

  SUBROUTINE get_SI_Ni_gradNk_triag( nfe, i, k, area, linearCoef, integral_out )
    ! Integral: Ni * grad(Nk) dA, OR Nk * grad(Ni) dA, over a triangle, with
    !     Ni as scalar shape function for the test point/node,
    !     Nk as the shape function for the second point/node (maybe the same as the test point if k==i)
    !     grad(Nk) : the gradient of Nk
    !
    ! Output: is a 2 by 1 vector
    !
    ! nfe: the polynomial degree (1: linear; 2: quadratic, 3: cubic etc)
    ! (Based on the formulae for high-order FE shape functions in Jin's book, 3rd edition)
    !
    ! linearCoef: the coefficients of 3 linear shape functions, pre-calculated
    !
    ! Jianbo Long, April, 2020
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nfe  ! degree of the polynomials (basis functions)
    INTEGER,INTENT(IN) :: i, k  ! local indices of shape functions within a triangle
    REAL(DPR),INTENT(IN) :: area, linearCoef(3,3)
    REAL(DPR),INTENT(OUT) :: integral_out(2)

    ! intermediate arguments
    REAL(DPR) :: t1, t2, t3, t4, grad1(2), grad2(2)
    INTEGER :: m, n, p, q
    TYPE(error_type) :: error
    CHARACTER(LEN=wordlen) :: msg
    CHARACTER(LEN=sub_len) :: thisSubroutine = "get_SI_Ni_gradNk_triag"

    SELECT CASE( nfe )
    CASE(1)
       ! the shape functions Ni,Nk are Linear shape functions
       IF( (i >= 1) .AND. (i <= 3) .AND. (k >= 1) .AND. (k <= 3) ) THEN
          CALL get_linear_SI_N123_triag( 1, 0, 0, area, t1 ) ! only Li(x,z)
          ! the gradient of Nk is contant
          integral_out = t1 *  linearCoef(k, 2:3)

       ELSE
          WRITE(msg,FMT='(A,A,A)') ' Invalid input of i and k (both should be >=1 and <=3) when <nfe = 1>: ',i,k
          CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
          IF(error_check(error)) CALL error_report(error)
       END IF

    CASE(2)       
       ! The shape functions are Quadratic order shape functions
       ! There are 4 situations regarding i & k
       IF( (i >= 1) .AND. (i <= 3) .AND. (k >= 1) .AND. (k <= 3) ) THEN
          ! situation 1: both points are vertex nodes
          grad1(1:2) = linearCoef(k, 2:3)

          IF( i /= k ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t1 ) ! Li^2 * Lk
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t3 ) ! Li * Lk
          ELSE
             CALL get_incomplete_linear_SI_N123_triag( 3, 0, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t3 )
          ENDIF

          CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t2 ) ! Li^2
          CALL get_incomplete_linear_SI_N123_triag( 1, 0, 0, t4 ) ! Li

          integral_out = (8 * t1 - 2 * t2 - 4 * t3 + t4) * (2.0 * area) * grad1

       ELSEIF( (i >= 1) .AND. (i <= 3) .AND. (k >= 4) .AND. (k <= 6) ) THEN
          ! situation 2: test node is a vertex, the other is a mid-edge node
          CALL get_two_vertices_from_mid_edge_point(k, m, n)
          grad1 = linearCoef(m, 2:3)  ! grad Lm
          grad2 = linearCoef(n, 2:3)

          ! m /= n, but one of (m,n) can be equal to i, or none of them is equal to i
          IF( (m == i) ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t1 ) ! Li^2 * Ln
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t2 ) ! Li *Ln
             CALL get_incomplete_linear_SI_N123_triag( 3, 0, 0, t3 ) ! Li^2 *Lm
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t4 ) ! Li *Lm
          ELSEIF( n == i ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 3, 0, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 0, 0, t2 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t3 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t4 )
          ELSEIF( (m/=i) .AND. (n/=i) )  THEN
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t1 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t2 )
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t3 )
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t4 )
          ENDIF

          integral_out = 4.0 * ( (2 * t1 - t2) * grad1 + (2 *t3 - t4) * grad2 ) * (2.0 * area)

       ELSEIF( (i >= 4) .AND. (i <= 6) .AND. (k >= 1) .AND. (k <= 3) ) THEN
          ! situation 3: test node is a mid-edge position, the other is a vertex
          CALL get_two_vertices_from_mid_edge_point(i, m, n)
          grad1 = linearCoef(k, 2:3)

          IF( (m == k) .OR. (n==k) ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t1 ) ! Lk*Lm*Ln
          ELSEIF( (m/=k) .AND. (n/=k) ) THEN
             CALL get_incomplete_linear_SI_N123_triag( 1, 1, 1, t1 )            
          ENDIF

          ! since m/=n anyway
          CALL get_incomplete_linear_SI_N123_triag( 1, 1, 0, t2 ) ! Lm*Ln
          integral_out = 4.0 * (4 * t1 - t2) * grad1 * (2.0 * area)

       ELSEIF( (i >= 4) .AND. (i <= 6) .AND. (k >= 4) .AND. (k <= 6) ) THEN
          ! situation 4: both nodes are mid-edge node
          ! The two edges can be : (1) the same; or (2) connected; in a triangle
          CALL get_two_vertices_from_mid_edge_point(i, m, n)  ! Ni(x,z), m/=n
          CALL get_two_vertices_from_mid_edge_point(k, p, q)  ! Nk(x,z), p/=q

          grad1 = linearCoef(p, 2:3)
          grad2 = linearCoef(q, 2:3)

          nes:IF( i == k ) THEN
             ! two edges are the same: m == p, and n == q
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t1 ) ! Lm*Ln*Lq
             CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t2 ) ! Lm*Ln*Lp
          ELSE
             ! two edges are connected, creating 4 scenarios further
             ! remember: m/=n, p/=q
             IF( (m == p) .OR. (n == p)) THEN
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 1, t1 ) ! Lm*Ln*Lq
                CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t2 ) ! Lm*Ln*Lp
             ELSEIF( (m == q) .OR. (n == q) ) THEN
                CALL get_incomplete_linear_SI_N123_triag( 2, 1, 0, t1 )
                CALL get_incomplete_linear_SI_N123_triag( 1, 1, 1, t2 )                   
             ENDIF
          ENDIF nes

          integral_out = 16.0 * (t1 * grad1 + t2 * grad2) * (2.0 * area)
       end IF

    CASE DEFAULT
       WRITE(msg,FMT='(A,A)') ' Invalid input of <nfe> (should be >=1 and <=2) by far: ', nfe
       CALL error_construct(error,ERROR_GENERAL,TRIM(ADJUSTL(thisModule)),TRIM(ADJUSTL(thisSubroutine)),msg)
       IF(error_check(error)) CALL error_report(error)
    end SELECT
    RETURN
  end SUBROUTINE get_SI_Ni_gradNk_triag

  ! -------------------------------------------------------------------
  ! -------------------------------------------------------------------

  SUBROUTINE get_linear_SI_N123_triag( norder1, norder2, norder3, area, integral_out )
    ! Integral: N1^n1 * N2^n2 * N3^n3 dV ,  with N1,N2,N3 as the 3 LINEAR shape functions
    !    within each triangle (associated with the 3 vertices, respectively);
    !    n1,n2,n3 as the power orders of the N1,N2,N3, respectively.
    !
    ! Based on the formula on eq(4.32), P85 in Jin's book, 3rd edition
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: norder1, norder2, norder3  ! power order of shape function
    REAL(DPR),INTENT(IN) :: area
    REAL(DPR),INTENT(OUT) :: integral_out

    CALL get_incomplete_linear_SI_N123_triag( norder1, norder2, norder3, integral_out )
    integral_out = integral_out * 2.0 * area
    RETURN
  end SUBROUTINE get_linear_SI_N123_triag
  !---------------------------------------------------------------
  !---------------------------------------------------------------

  SUBROUTINE get_incomplete_linear_SI_N123_triag( norder1, norder2, norder3, integrand )
    ! Integral: N1^n1 * N2^n2 * N3^n3 dV ,  with N1,N2,N3 as the 3 LINEAR shape functions
    !    within each triangle (associated with the 3 vertices, respectively);
    !    n1,n2,n3 as the power orders of the N1,N2,N3, respectively.
    !
    ! Based on the formula on eq(4.32), P85 in Jin's book, 3rd edition
    !
    ! ONLY get the integration results without the coefficient '2*A', A is the area of the triangle.
    !   (hence, the 'incomplete' integration; the purpose is to reduce float point operations)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: norder1, norder2, norder3  ! power order of shape function
    REAL(DPR),INTENT(INOUT) :: integrand
    INTEGER :: numerator, denor

    numerator = get_factorial( norder1 ) * get_factorial( norder2 ) * get_factorial( norder3 )

    denor = get_factorial( 2 + norder1 + norder2  + norder3 )
    ! convert both operators to float numbers before division (otherwise can cause modulo operations)
    integrand = (numerator * 1.0) / (denor * 1.0)
    RETURN
  end SUBROUTINE Get_Incomplete_Linear_SI_N123_Triag

  !---------------------------------------------------------------
  !---------------------------------------------------------------

  SUBROUTINE get_basisWeights_linear_triag( nnode, linearCoef, area)
    ! For 2D triangular elements
    ! --get the expansion coefficients(i.e., a,b,c in a + bx + cz) for
    ! 3 linear basis functions N1,N2, and N3. This process is essentially
    ! asking for the inverse of the 3 by 3 linear interpolation matrix.
    !
    ! -- Once the basis functions are available, then an approximant within the element
    !   can be written as: f(x,z) = N1*f1 + N2*f2 + N3 * f3
    !
    ! on output: area and the coefficients of shape functions
    !               C1   C2  C3 (Columns)   
    ! N1:           a    b    c
    ! N2:           -    -    -
    ! N3:           -    -    -
    !(rows)
    ! Based on formulae on P84 in Jin's book, 3rd edition.

    USE derived_data_module, ONLY: nodes
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nnode(3)   ! 3 global nodal numbers of the element. The order
    !  of the numbering is taken as that from the mesh files.
    REAL(DPR),INTENT(INOUT) :: linearCoef(3, 3), area
    REAL(DPR) :: x1, x2, x3, z1, z2, z3
    REAL(DPR) :: a1,b1,c1, a2,b2,c2, a3,b3,c3, p1(2), p2(2), p3(2)

    x1 = nodes( nnode(1) )%x;  z1 = nodes( nnode(1) )%z

    x2 = nodes( nnode(2) )%x;  z2 = nodes( nnode(2) )%z

    x3 = nodes( nnode(3) )%x;  z3 = nodes( nnode(3) )%z

    ! vector form of the positions
    p1 = [x1, z1];  p2 = [x2, z2];  p3 = [x3, z3]

    a1 = x2 * z3 - z2 * x3;  b1 = z2 - z3;  c1 = x3 - x2
    a2 = x3 * z1 - z3 * x1;  b2 = z3 - z1;  c2 = x1 - x3
    a3 = x1 * z2 - z1 * x2;  b3 = z1 - z2;  c3 = x2 - x1

    CALL get_area_general_triag(p1, p2, p3, area)

    ! N_i(x,y,z) = a_i + b_i * x + c_i * z;  i = 1, 2, 3

    ! all "a"s in 3 basis functions
    linearcoef(:, 1) = [a1, a2, a3] / (2.0 * area)

    ! all "b"s
    linearcoef(:, 2) = [b1, b2, b3] / (2.0 * area)

    ! all "c"s
    linearcoef(:, 3) = [c1, c2, c3] / (2.0 * area)
    RETURN
  end SUBROUTINE get_basisWeights_linear_triag

  SUBROUTINE get_shapeFunc_value_triag(nshape, BasisWeights, pos, fun, funDx, funDz)
    USE modelling_parameter, ONLY: modelling
    IMPLICIT NONE
    ! shape functions' (and their derivatives') values at a particular 2D position; this is to facilitate
    ! FE-based post-processing
    INTEGER,INTENT(IN) :: nshape   ! # of shape functions
    INTEGER :: k
    REAL(DPR),INTENT(IN) :: BasisWeights(3,3), pos(2)   ! ceofs of linear shape funcs; 2D location
    REAL(DPR),INTENT(INOUT) :: fun(nshape), funDx(nshape), funDz(nshape)
    REAL(DPR) :: L(3), g1(2), g2(2), g3(2), gf(2)

    IF(modelling%FEdegree == 1) THEN
       DO k = 1, nshape
          ! a_k + b_k*x0 + c_k*z0 = N_k(x0,z0)
          fun(k) = BasisWeights(k, 1) + BasisWeights(k, 2) * pos(1) +  BasisWeights(k, 3) * pos(2)
          funDx(k) = BasisWeights(k, 2)   ! partial derivative: df/dx
          funDz(k) = BasisWeights(k, 3)
       end DO
    ELSEIF(modelling%FEdegree == 2) THEN
       ! nshape = 6
       ! a_k + b_k*x0 + c_k*z0 = N_k(x0,z0)
       DO k = 1, 3
          L(k) = BasisWeights(k, 1) + BasisWeights(k, 2) * pos(1) +  BasisWeights(k, 3) * pos(2)
          fun(k) = 2.0 * L(k)**2 - L(k)
          funDx(k) = (4.0 * L(k) - 1.0) * BasisWeights(k, 2)
          funDz(k) = (4.0 * L(k) - 1.0) * BasisWeights(k, 3)
       END DO

       g1 = [BasisWeights(1, 2), BasisWeights(1, 3)]  ! grad L1
       g2 = [BasisWeights(2, 2), BasisWeights(2, 3)]  ! grad L2
       g3 = [BasisWeights(3, 2), BasisWeights(3, 3)]  ! grad L3
       ! 4th shape func
       fun(4) = 4.0 * L(1) * L(2)
       gf = 4.0 * ( g1 * L(2) + L(1)* g2 )
       funDx(4) = gf(1)
       funDz(4) = gf(2)
       ! 5th shape func
       fun(5) = 4.0 * L(2) * L(3)
       gf = 4.0 * ( g2 * L(3) + L(2)* g3 )
       funDx(5) = gf(1)
       funDz(5) = gf(2)

       ! 6th shape func
       fun(6) = 4.0 * L(1) * L(3)
       gf = 4.0 * ( g1 * L(3) + L(1)* g3 )
       funDx(6) = gf(1)
       funDz(6) = gf(2)
    end IF

  end SUBROUTINE get_shapeFunc_value_triag

  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------
  SUBROUTINE get_area_general_triag(p1, p2, p3, tri_area)
    ! get the area of a general triangle in 2-D based on its 3 vertices positions
    !  for point i (i=1,2,3):
    !  pi(1) -- x coordinate
    !  pi(2) -- z ..
    ! area = half of the absolute value of the determinant of the 3 by 3 matrix
    IMPLICIT NONE

    REAL(DPR),INTENT(IN)  :: p1(2), p2(2), p3(2)
    REAL(DPR),INTENT(OUT) :: tri_area     ! surface area of the element

    REAL(DPR)  :: b1, c2, b2, c1

    b1 = p2(2) - p3(2)   ! z2 - z3
    c2 = p1(1) - p3(1)   ! x1 - x3
    b2 = p3(2) - p1(2)
    c1 = p3(1) - p2(1)   ! x3 - x2

    tri_area = (1.0 / 2.0) * ABS( b1 * c2 - b2 * c1 )

    RETURN
  end SUBROUTINE get_area_general_triag
  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------
  SUBROUTINE get_barycentric_coord_general_triag(pt, p1, p2, p3, bary, errFlag )
    ! get the barycentric coordinates of an arbitrary 2D point based on the 3 vertex positions
    !  of the triangle under consideration; the coordinates of the 3 vertices of the
    !  triangle are given as input:
    !  for i=1,2,3
    !  pi(1) -- x coordinate of the i-th vertex
    !  pi(2) -- z
    !
    !  pt(1:2): the Cartesian coordinates of the test point in 2D
    ! -- On output:
    !  bary(3): the barycentric coordinates (l1, l2, l3) of the test point
    !    (with L1+L2 +L3 = 1; 0<Li<1,i=1,2,3 if the test point lies inside )
    !  errFlag: error flag of the calculation
    !    = 0 : the triangle is not degenerated
    !    > 0 : the input triangle is degenerated (e.g., a point, or a segment) in a numerical sense
    ! -- Based on the discussion in Jin's book, third edn, (4.169-175), P133.
    !  The barycentric coordinates, or areal coordinates of a point are ESSENTIALLY the values
    !    of 3 linear shape functions at the test point.
    IMPLICIT NONE

    REAL(DPR),INTENT(IN)  :: p1(2), p2(2), p3(2), pt(2)
    REAL(DPR),INTENT(OUT) :: bary(3)
    INTEGER,INTENT(OUT) :: errFlag

    REAL(DPR) :: area
    REAL(DPR) :: x1, x2, x3, z1, z2, z3, x, z
    REAL(DPR) :: a1,b1,c1, a2,b2,c2, a3,b3,c3
    REAL(DPR) :: zeroTol = 1.d-40

    x1 = p1(1);  z1 = p1(2)
    x2 = p2(1);  z2 = p2(2)
    x3 = p3(1);  z3 = p3(2)

    x = pt(1);  z = pt(2)

    a1 = x2 * z3 - z2 * x3;  b1 = z2 - z3;  c1 = x3 - x2
    a2 = x3 * z1 - z3 * x1;  b2 = z3 - z1;  c2 = x1 - x3
    a3 = x1 * z2 - z1 * x2;  b3 = z1 - z2;  c3 = x2 - x1

    CALL get_area_general_triag(p1, p2, p3, area) ! triangle's area

    IF( ABS(area) <= zeroTol )  THEN
       errFlag = 1  ! degenerated triangle
       bary = 0.0
       RETURN
    ELSE
       bary(1) = (a1 + b1 * x + c1 * z) / (2.0 * area)
       bary(2) = (a2 + b2 * x + c2 * z) / (2.0 * area)
       bary(3) = (a3 + b3 * x + c3 * z) / (2.0 * area)
       errFlag = 0
    END IF
    RETURN
  end SUBROUTINE get_barycentric_coord_general_triag
  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------

  SUBROUTINE whether_point_inside_triag(pt, p1, p2, p3, ins, errFlag  )
    ! test whether a point is inside a triangle or not; the coordinates of the 3 vertices of the
    !  triangle are given as input:
    !  for i=1,2,3
    !  pi(1) -- x coordinate of the i-th vertex
    !  pi(2) -- z
    !
    !  pt: the Cartesian coordinates of the test point
    ! -- On output:
    !  ins: the integer indicator of the situation
    !    = -1: an error occurs
    !    = 0 : strictly outside
    !    = 1 : point at one of the vertices
    !    = 2 : point on one of the facets (edges) of the triangle
    !    = 3 : point strictly inside the triangle

    !  errFlag: error flag of the calculation
    !    = 0 : the triangle is not degenerated
    !    > 0 : the input triangle is degenerated (e.g., a point, or a segment)
    IMPLICIT NONE
    REAL(DPR),INTENT(IN)  :: p1(2), p2(2), p3(2), pt(2)
    INTEGER,INTENT(INOUT) :: ins, errFlag
    INTEGER :: k
    REAL(DPR) :: bary(3), x(3), z(3)

    ! used to detect if two points are the same in 2-D space
    REAL(DPR) :: zeroTol = 1.d-3

    x = [p1(1), p2(1), p3(1)]
    z = [p1(2), p2(2), p3(2)]

    ! check if the input point is one of the vertices
    DO k = 1, 3
       IF( ABS(pt(1)-x(k)) <= zeroTol .AND. ABS(pt(2)-z(k)) <= zeroTol ) THEN
          ins = 1
          errFlag = 0
          RETURN
       end IF
    end DO

    ! get the barycentric coordinates
    CALL get_barycentric_coord_general_triag(pt, p1, p2, p3, bary, errFlag )

    IF( errFlag > 0 )  THEN
       ! element degenerated
       ins = -1
       RETURN
    ELSE
       !IF( MINVAL(bary) >= 0.0  .AND. MAXVAL(bary) <= 1.0 ) THEN
       ! small value here to compensate for the numerical error in the case of exactly
       ! zero areal coordinate of degenerated sub-triangle
       IF( MINVAL(bary) >= -1.d-3  .AND. MAXVAL(bary) <= 1.0 ) THEN  
          ! when there is a zero barycentric coordinate
          IF( ANY( ABS(bary) <= 1.D-20 ) )  THEN
             ins = 2
          ELSE
             ins = 3
          ENDIF
       ELSE
          ins = 0
       ENDIF
    END IF

    RETURN
  end SUBROUTINE Whether_Point_Inside_Triag

  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------

  SUBROUTINE which_triag_for_point( pt, iele, errFlag )
    ! the tree structure of the point data is used to find the nearest node of the mesh
    !   for a given point
    USE derived_data_module, ONLY: nodes
    USE mesh2d_discretization, ONLY: ele2node, p2ele
    IMPLICIT NONE

    REAL(DPR),INTENT(IN) :: pt(2)   ! the inquiry point
    INTEGER,INTENT(OUT) :: iele  ! the # of the found element in the mesh
    INTEGER,INTENT(OUT) :: errFlag  ! error flag; =0:No errors; >0: error happens

    INTEGER :: i, j, nele, test_ele, Points(3)
    INTEGER :: inside, error
    REAL(DPR)  :: p1(2), p2(2), p3(2)
    REAL(DPR) :: x(3), z(3)
    LOGICAL :: isFound

    CHARACTER(LEN=sub_len) :: thisSub = "which_triag_for_point"
    REAL(DPR) :: rs, R0

    isFound = .FALSE.
    errFlag = 0

    IF( (isFound .EQV. .FALSE.) ) THEN
       ! switch to the global search (very slow)...
       !PRINT*, '------DEBUG TESTING ---------'
       !PRINT*,''
       !PRINT*,'......debugging....USED global search for the point in a tetrahedron ...'
       nele = SIZE( ele2node, 1 )
       DO i = 1, nele
          test_ele = i
          Points = ele2node( test_ele, 1 : 3 )
          DO j = 1, SIZE(Points)
             x(j) = nodes( Points(j) )%x
             z(j) = nodes( Points(j) )%z
          end DO
          p1 = [x(1), z(1)]
          p2 = [x(2), z(2)]
          p3 = [x(3), z(3)]
          CALL whether_point_inside_triag(pt, p1, p2, p3, inside, error  )

          IF( inside >= 1 ) THEN
             ! the element is found
             iele = test_ele
             isFound = .TRUE.
             RETURN
          END IF
       END DO
    end IF

    IF( (isFound .EQV. .FALSE.) ) THEN
       PRINT*, '----- Warning Error ------'
       PRINT*, 'In subroutine <'//TRIM(ADJUSTL(thisSub))//'>'
       PRINT*, 'Fatal error in searching the element for a given point position:'
       PRINT*, 'No element was found for the input point: ', pt
       PRINT*,'The element containing this point is possibly degenerated within the given numerical accuracy!'
       PRINT*, '-------ending warning-----'
       ! choose the element with the closest node in the mesh for the point as the one element
       DO i = 1, SIZE(nodes)
          p1 = [nodes(i)%x, nodes(i)%z]
          rs = SQRT( SUM((p1 - pt)**2) )
          IF (i == 1) THEN
             R0 = rs
             CYCLE
          ELSE
             IF( rs <= R0 ) THEN
                R0 = rs
                iele = p2ele(i)%ele(1)
                !errFlag = 1
             end IF
          END IF
       END DO
    end if
    RETURN
  end SUBROUTINE which_triag_for_point

  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------

  PURE FUNCTION get_factorial( kn ) RESULT(f)
    ! get the factorial k! for non-negative integer
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: kn
    INTEGER :: j
    INTEGER  :: f
    f = 1
    IF( kn > 0 ) THEN
       DO j = 1, kn
          f = f * j
       END DO
    ELSE
       f = 1
    END IF
    RETURN
  end FUNCTION Get_Factorial

  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------
  SUBROUTINE find_map_array1_to_array2( n1, index1, n2, index2, mapIndex )
    ! Purpose: to scatter elements of array1 into corresponding positions in array2.
    !  For example, i-th element in array1 goes to j-th element in array2 if
    !    index1(i) = index2(j)
    !
    ! E.g.,  index1 entries: /4, 3,  10,  7/
    !        index2 entries: /4, 22, 9, 10, 7, 6, 3, 12/ (contains all entries in index1)
    !
    ! then the output map index from 1 to 2 should be: (1, 7, 4, 5)
    ! i.e., entry 3 in index1 is the second element in index1, but the 7-th element in index2.

    ! This is used to improve efficiency in assembling coeffs from all connected elements
    implicit none

    INTEGER,INTENT(IN) :: n1, n2, index1(n1),  index2(n2)
    INTEGER,INTENT(OUT) :: mapIndex( n1 )  ! stores index values in the second array
    !                                         where its entries will be filled by the
    !                                         first array entries.
    INTEGER :: k, j

    mapIndex = 0
    DO k = 1, n1
       innerLoop: DO j = 1, n2
          IF( index1(k) == index2(j) )  THEN
             mapIndex( k ) = j
             EXIT innerLoop
          END IF
       END DO innerLoop
    END DO
    RETURN
  end SUBROUTINE find_map_array1_to_array2

  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------

  SUBROUTINE mapping_array1_to_array2( n1, array1, mapIndex, n2, array2 )
    ! Purpose: to scatter elements of array1 into corresponding positions in array2.
    !  For example, i-th element in array1 goes to j-th element in array2 if
    !    index1(i) = index2(j)
    !
    ! This is used to assemble coefficients from all connected elements
    implicit none
    INTEGER,INTENT(IN) :: n1, n2, mapIndex(n1)
    REAL(DPR),INTENT(IN) :: array1(n1)
    REAL(DPR),INTENT(OUT) :: array2(n2)
    INTEGER :: k
    array2 = 0.d0
    DO k = 1, n1
       array2( mapIndex(k) ) = array1( k )
    END DO
    RETURN
  end SUBROUTINE mapping_array1_to_array2

  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------

  SUBROUTINE get_two_vertices_from_mid_edge_point(k, m, n)
    ! --For triangles
    !
    ! For high-order FE methods, get hierarchical two vertex points based on the input mid-edge point
    ! The numebering of non-vertex nodes follow standard hierarchical order, e.g., in Jin's book
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: k  ! the local index of the mid-edge point
    INTEGER :: m, n

    m = 0;  n = 0

    SELECT CASE (k)
    CASE(4)
       m = 1;  n = 2
    CASE(5)
       m = 2;  n = 3
    CASE(6)
       m = 1;  n = 3
    END SELECT
    RETURN
  end SUBROUTINE get_two_vertices_from_mid_edge_point
  !---------------------------------------------------------------
  !---------------------------------------------------------------
  SUBROUTINE connected_edges_for_vert_search( vert, list )
    ! To find all the connecting edge-type dofs (mid-edge points) for
    !   a test vertex; these dofs all belong to the cells connecting the test vertex.
    ! --required for high-order FE
    ! -- triangular cells type
    ! Output:
    !  list: a list of all connected edge-type dofs, WITHOUT the test vertex
    USE mesh2d_discretization, ONLY: p2ele, ele2edge
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: vert  ! the (global index of) test vertex
    INTEGER,INTENT(INOUT),ALLOCATABLE ::   list(:)
    INTEGER   :: nele, k, j, m, currentEdge, currentEle

    nele = SIZE(p2ele(vert)%ele)  ! number of connected cells for the test vertex
    m = 0
    DO k = 1, nele
       currentEle = p2ele(vert)%ele(k)
       DO j = 1, 3     ! each cell (triangle) has 3 edges
          m = m + 1
          currentEdge = ele2edge( currentEle, j)
          IF( m == 1 ) THEN
             CALL add_to_list_integer( list, currentEdge)
          ELSE
             IF( .NOT. ANY( list == currentEdge) ) CALL add_to_list_integer( list, currentEdge)
          ENDIF
       ENDDO
    ENDDO

    RETURN
  end SUBROUTINE connected_edges_for_vert_search

  SUBROUTINE connected_verts_for_edge_search( iedge, list )
    ! To find all the connecting vertex-type dofs for
    !   a test edge; these dofs all belong to the cells connecting the test edge.
    ! --required for high-order FE
    ! -- triangular cells type
    ! Output:
    !  list: a list of all connected vertex-type dofs, WITHOUT the test edge
    USE mesh2d_discretization, ONLY: eg2ele, ele2node
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: iedge  ! the (global index of) test edge
    INTEGER,INTENT(INOUT),ALLOCATABLE ::   list(:)
    INTEGER   :: nele, k, j, m, currentVert, currentEle

    nele = SIZE(eg2ele(iedge)%ele)  ! number of connected cells for the test edge
    m = 0
    DO k = 1, nele
       currentEle = eg2ele(iedge)%ele(k)
       DO j = 1, 3     ! each cell (triangle) has 3 vertices
          m = m + 1
          currentVert = ele2node( currentEle, j)
          IF( m == 1 ) THEN
             CALL add_to_list_integer( list, currentVert)
          ELSE
             IF( .NOT. ANY( list == currentVert) ) CALL add_to_list_integer( list, currentVert)
          ENDIF
       ENDDO
    ENDDO

    RETURN
  end SUBROUTINE connected_verts_for_edge_search

end MODULE scalar_FE_kernels_2d
