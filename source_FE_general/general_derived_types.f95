
MODULE float_precision
  IMPLICIT NONE
  ! to make the range and precision of float point numbers more portable
  INTEGER,PARAMETER :: DPR = SELECTED_REAL_KIND(15,307)  ! i.e., to keep 15 sig digits, and the range is 10^{-307} to 10^{307}.
  INTEGER,PARAMETER :: CDPR = SELECTED_REAL_KIND(15,307)  ! long complex type
  INTEGER,PARAMETER :: SPR = SELECTED_REAL_KIND(7,37)   ! the usual single precision real value
  
END MODULE float_precision

MODULE derived_data_module

  USE float_precision, ONLY: DPR, CDPR
  IMPLICIT NONE
  ! Declare type node_info;

  TYPE:: geometry_node_type
     ! basic mesh-related geometric node type
     INTEGER    :: bud     ! boundary info: 
     !  bud == 1 : on the computational outer boundaries;
     !  bud /= 1 : inside computational outer boundaries, in this case,
     !             value of 'reg' equals to the attribute value of the region
     !             in which the node is located.
     INTEGER,ALLOCATABLE    :: reg(:)       ! regional attribution info
     REAL(DPR)  :: x       ! x-coordinate
     REAL(DPR)  :: y       ! y-coordinate
     REAL(DPR)  :: z       ! z-coordinate
  end type geometry_node_type

  TYPE,EXTENDS(geometry_node_type):: geoph_grav_node_type
     REAL(DPR)  :: rho     ! physical property, density (kg/m^3)
     REAL(DPR)  :: phig     ! gravity potential data (m^2/ s^2), or other real-valued scalar potential (e.g., DC)
  end TYPE geoph_grav_node_type


  TYPE,EXTENDS(geoph_grav_node_type) :: geoph_EM_node_type
     REAL(DPR)  :: sigma     ! physical property, conductivity
     REAL(DPR)  :: mu        ! physical property, magnetic permeability
     REAL(DPR)  :: permt     ! physical property, electrical permittivity

     COMPLEX(CDPR)  :: A(3)  !  vector potential value associated
     COMPLEX(CDPR)  :: phi   ! EM scalar potential
     COMPLEX(CDPR)  :: E(3)   ! scalar components of Electric field
     COMPLEX(CDPR)  :: H(3)   ! scalar components of magnetic field
  END type Geoph_EM_Node_Type

  TYPE(geometry_node_type),ALLOCATABLE :: nodes(:)
  TYPE(geoph_EM_node_type),ALLOCATABLE :: dofs(:)

CONTAINS
  SUBROUTINE deep_node_copy(nod1, nod2)
    ! COPY nod1 TO nod2
    IMPLICIT NONE
    TYPE(geoph_EM_node_type),INTENT(IN) :: nod1
    TYPE(geoph_EM_node_type),INTENT(INOUT) :: nod2
    nod2%x = nod1%x
    nod2%y = nod1%y
    nod2%z = nod1%z
    nod2%reg = nod1%reg
    nod2%bud = nod1%bud
    nod2%rho = nod1%rho
    nod2%sigma = nod1%sigma
    nod2%mu = nod1%mu
    nod2%permt = nod1%permt
    nod2%phig = nod1%phig
    nod2%A = nod1%A
    nod2%phi = nod1%phi
    nod2%E = nod1%E
    nod2%H = nod1%H
    RETURN
  end SUBROUTINE deep_node_copy

  SUBROUTINE deep_node_copy_geometry2EM(nod1, nod2)
    ! COPY nod1 of geometry type TO nod2 of EM type
    IMPLICIT NONE
    TYPE(geometry_node_type),INTENT(IN) :: nod1
    TYPE(geoph_EM_node_type),INTENT(INOUT) :: nod2
    nod2%x = nod1%x
    nod2%y = nod1%y
    nod2%z = nod1%z
    nod2%reg = nod1%reg
    nod2%bud = nod1%bud
    RETURN
  end SUBROUTINE deep_node_copy_geometry2EM


END MODULE Derived_Data_Module



MODULE constants_module
  ! 
  ! Purpose: assign physical or mathematical constants in a finite-precison type
  !
  USE float_precision, ONLY: DPR
  IMPLICIT NONE

  ! Newton's gravitational constant, SI unit: m^3.Kg^-1.s^-2
  REAL(DPR),PARAMETER  :: gravity_cst = 6.67408d-11

  REAL(DPR),PARAMETER  :: pi = 3.14159d0

  REAL(DPR),PARAMETER  :: mu0 = 4 * pi * 1.D-7

  REAL(DPR),PARAMETER  :: epsilon0 = 8.854187817 * 1.D-12


  COMPLEX(DPR),PARAMETER :: cmpx_one = CMPLX( 1.d0, 0.d0, KIND= DPR )

  COMPLEX(DPR),PARAMETER :: cmpx_i = CMPLX( 0.d0, 1.d0, KIND= DPR )

  COMPLEX(DPR),PARAMETER :: cmpx_zero = CMPLX( 0.d0, 0.d0, KIND= DPR )

  ! unit normal vector in 3D Cartesian coordinates
  REAL(DPR),PROTECTED :: unit_vec_x(3) = [ 1.d0, 0.d0, 0.d0 ]
  REAL(DPR),PROTECTED :: unit_vec_y(3) = [ 0.d0, 1.d0, 0.d0 ]
  REAL(DPR),PROTECTED :: unit_vec_z(3) = [ 0.d0, 0.d0, 1.d0 ]
  
END MODULE constants_module


