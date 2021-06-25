MODULE EM_RHS_concatenate_module
  ! Author: Jianbo Long, Mar, 2020
  USE mesh_concatenate_module, ONLY: ndof
  USE modelling_parameter, ONLY: modelling
  IMPLICIT NONE
  PRIVATE
  INTEGER,PARAMETER :: sub_len = 300
  !CHARACTER(LEN=sub_len),PROTECTED :: thisModule = "EM_RHS_concatenate_module"
  PUBLIC :: EM_RHS_concatenate

CONTAINS
  ! -----------------------------------------------------------------------

  SUBROUTINE EM_RHS_concatenate( )
    ! Purposes: Overall interface procedure .
    !USE float_precision, ONLY: CDPR
    USE constants_module, ONLY: cmpx_zero !, cmpx_i
    USE FORTRAN_generic,ONLY: print_general_msg
    USE linear_system_equation_data, ONLY: initialize_matrices_EM_RHS, EM_RHS, EM_solution, &
         initialize_matrices_gt_RHS, gt_RHS, gt_solution
    IMPLICIT NONE

    IF( TRIM(ADJUSTL(modelling%DataType)) == 'CM' )  THEN
       CALL initialize_matrices_EM_RHS()
       CALL print_general_msg('Complex-data Right-hand-sides')
       CALL print_general_msg('PDE: Helmholtz for Ey and Hy')
       ALLOCATE( EM_RHS(ndof, 2) )    ! (Ey, Hy)
       ALLOCATE( EM_solution(ndof, 2) )
       EM_RHS = cmpx_zero
       CALL get_boundary_from_files_DD_test(EM_RHS(:, 1), EM_RHS(:, 2))
    ELSE
       CALL initialize_matrices_gt_RHS()
       CALL print_general_msg('Real-data Right-hand-sides')
       ALLOCATE( gt_RHS(ndof, 1) )
       ALLOCATE( gt_solution(ndof, 1) )
       gt_RHS = 0.0
       CALL get_boundary_from_files_gt_data(gt_RHS(:,1))
    END IF

    RETURN
  END SUBROUTINE EM_RHS_Concatenate

  ! -----------------------------------------------------------------------
  SUBROUTINE get_boundary_from_files_gt_data(RHS)
    USE float_precision, ONLY: DPR
    USE FORTRAN_generic,ONLY: print_general_msg
    IMPLICIT NONE
    REAL(DPR),INTENT(INOUT) :: RHS(:)
    REAL(DPR),ALLOCATABLE :: rad(:)
    !REAL(DPR),ALLOCATABLE :: x(:), z(:)
    INTEGER :: k, null, ndata, fid

    CALL print_general_msg('read files of boundary values')

    fid = 48

    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//'All_rad.txt', STATUS = 'OLD')
    READ(fid, *) ndata
    ALLOCATE(rad(ndata))
    DO k = 1, ndata
       READ(fid, *) null, rad(k)
    end DO
    CLOSE(UNIT=fid)



    RETURN
  end SUBROUTINE get_boundary_from_files_gt_data
  ! -----------------------------------------------------------------------
  SUBROUTINE get_boundary_from_files_DD_test(RHS_Ey, RHS_Hy)
    USE float_precision, ONLY: CDPR, DPR
    USE derived_data_module, ONLY: nodes
    USE file_rw, ONLY: read_sol_Ron_MT2D
    USE FORTRAN_generic,ONLY: print_general_msg
    IMPLICIT NONE
    COMPLEX(CDPR),INTENT(INOUT) :: RHS_Ey(:), RHS_Hy(:)
    COMPLEX(CDPR),ALLOCATABLE :: solEy(:), solHy(:) ! global solutions
    REAL(DPR),ALLOCATABLE :: x(:), z(:)
    INTEGER,ALLOCATABLE :: nodemap(:)
    INTEGER :: k 

    CALL print_general_msg('Ron-DD: read files of Ey and Hy')
    CALL read_sol_Ron_MT2D(solEy, solHy, x, z, nodemap )
    ! for debugging
    IF( SIZE(RHS_Ey) /= SIZE(nodemap) ) THEN
       PRINT*, 'nodemap AND RHS_Ex not matching !!'
       PRINT*, "SIZE(RHS_Ey) = ", SIZE(RHS_Ey)
       PRINT*, "SIZE(nodemap) = ", SIZE(nodemap)
       STOP
    END IF

    DO k = 1, SIZE(RHS_Ey)
       IF(nodes(k)%bud == 1) THEN
          RHS_Ey(k) = solEy( nodemap(k) )
          RHS_Hy(k) = solHy( nodemap(k) )
       END IF
    END DO
    RETURN
  end SUBROUTINE get_boundary_from_files_DD_test


end MODULE EM_RHS_Concatenate_Module






