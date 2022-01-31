PROGRAM mesh2d_partition
  ! for partitioning a 2-D mesh into submeshes.
  ! Author: Jianbo Long, April, 2020
  USE float_precision, ONLY: DPR
  USE FORTRAN_generic, ONLY: percentage_display
  USE FORTRAN_generic, ONLY: add_to_list_integer, find_intersect_integer, print_general_msg,&
       print_debug_msg, print_error_msg

  IMPLICIT NONE

  ! basic info for unstructured meshes
  INTEGER,ALLOCATABLE :: ele2node(:,:), neighele(:,:)  
  !!INTEGER,ALLOCATABLE :: edge2node(:,:)  !!, ele2edge(:,:)

  TYPE :: node_cls
     REAL(DPR)  :: x       ! x-coordinate
     !REAL(DPR)  :: y       ! y-coordinate
     REAL(DPR)  :: z       ! z-coordinate
     INTEGER    :: bud       ! boundary info
     INTEGER,ALLOCATABLE  :: reg(:)       ! regional attribution info
  end type node_cls

  TYPE(node_cls),ALLOCATABLE :: nodes(:), new_nodes(:)
  INTEGER,ALLOCATABLE :: nodemap(:), node_g2l_map(:),new_ele2node(:,:), new_neighele(:,:)
  REAL(DPR) :: tinyValue = 1.D-3
  INTEGER,PARAMETER :: nregion = 3
  INTEGER,PARAMETER :: regionMark(3) = (/10, 20, 30/)
  CHARACTER(LEN=200),PARAMETER :: basename = "MT_2D_3layer.1"
  CHARACTER(LEN=400),PARAMETER :: mesh_path = "/media/jack/NewVolume/3DEM_Jianbo_test_data/MT2D_meshes/"
  INTEGER :: k, mark


  CALL read_2D_mesh_from_files(nregion, regionMark, basename, mesh_path)

  DO k = 1, nregion
     mark = regionMark(k)
     CALL get_sub_nodes(mark, nodes, new_nodes, node_g2l_map, nodemap)
          PRINT*, 'K=', K
     CALL get_sub_elements(mark, ele2node, neighele, node_g2l_map, new_ele2node, new_neighele)

     CALL write_mesh_files(k, mesh_path, basename, new_nodes, new_ele2node, new_neighele, nodemap, node_g2l_map)
     
     IF(ALLOCATED(new_nodes)) DEALLOCATE(new_nodes)
     IF(ALLOCATED(nodemap)) DEALLOCATE(nodemap)
     IF(ALLOCATED(node_g2l_map)) DEALLOCATE(node_g2l_map)
     IF(ALLOCATED(new_ele2node)) DEALLOCATE(new_ele2node)
     IF(ALLOCATED(new_neighele)) DEALLOCATE(new_neighele)

  end DO


CONTAINS

  !! ------------------------------------------------------------------------------------------
  SUBROUTINE read_2D_mesh_from_files(nregion, regionMark, basename, mesh_path)
    USE FORTRAN_fileio_mod, ONLY : fileio_exists
    IMPLICIT NONE

    INTEGER,INTENT(IN) :: nregion, regionMark(:)
    CHARACTER(LEN=*),INTENT(IN) :: basename, mesh_path
    ! INput files
    CHARACTER(LEN=200)  :: nodefile  
    CHARACTER(LEN=200)  :: elefile  
    !CHARACTER(LEN=200)  :: edgefile 
    CHARACTER(LEN=200)  :: neighfile


    nodefile = TRIM(ADJUSTL(basename))//".node"
    elefile = TRIM(ADJUSTL(basename))//".ele"
    !edgefile = TRIM(ADJUSTL(basename))//".edge"
    neighfile = TRIM(ADJUSTL(basename))//".neigh"

    PRINT*,'--------------------------------------------------------------------------------'
    PRINT*,'-----------        mesh files found       --------------------------------------'
    PRINT*,'--------------------------------------------------------------------------------'
    PRINT*,TRIM(ADJUSTL(mesh_path))
    IF(fileio_exists( TRIM(ADJUSTL(mesh_path))//TRIM(ADJUSTL(nodefile)) )) PRINT*,TRIM(ADJUSTL(nodefile))
    IF(fileio_exists( TRIM(ADJUSTL(mesh_path))//TRIM(ADJUSTL(elefile)) )) PRINT*,TRIM(ADJUSTL(elefile))
    IF(fileio_exists( TRIM(ADJUSTL(mesh_path))//TRIM(ADJUSTL(neighfile)) )) PRINT*,TRIM(ADJUSTL(neighfile))

    CALL read_node_file(mesh_path, nodefile, nodes)
    CALL check_node_boundary(nodes)
    CALL read_element_file(mesh_path, elefile, regionMark, ele2node)
    CALL update_node_regional_marker(nregion, regionMark)

    CALL read_neigh_element_file(mesh_path, neighfile, neighele)    

    RETURN
  end SUBROUTINE read_2D_mesh_from_files
  ! ---------------------------------------

  SUBROUTINE read_node_file(filepath, filename, nodeobj)
    IMPLICIT NONE

    INTEGER :: fid, npoint, i, inull
    CHARACTER(LEN=*),INTENT(IN) :: filepath, filename
    TYPE(node_cls),INTENT(INOUT),ALLOCATABLE :: nodeobj(:)

    fid = 41
    ! read .node file
    OPEN(UNIT=fid, FILE= TRIM(ADJUSTL(filepath))//TRIM(ADJUSTL(filename)), STATUS = 'OLD')
    READ(fid,*) npoint

    IF(ALLOCATED(nodeobj) .EQV. .TRUE.)  DEALLOCATE( nodeobj );  ALLOCATE( nodeobj(npoint) )

    DO i = 1, npoint
       ! In 2-D Triangle files, boundary marker is 0, -1, or 1, or other user-defined values.
       ! In 2-D mesh, y is undefined.
       READ(fid,*)  inull, nodeobj(i)%x, nodeobj(i)%z, nodeobj(i)%bud
    END DO
    CLOSE(fid)
  end SUBROUTINE read_node_file

  ! ---------------------------------------

  SUBROUTINE check_node_boundary(nodeobj)
    ! here, "node" and "point" refer to the same thing.
    ! -- ONLY applicable for rectangular domains !
    IMPLICIT NONE

    INTEGER :: npoint, i
    TYPE(node_cls),INTENT(INOUT),ALLOCATABLE :: nodeobj(:)
    REAL(DPR) :: min_value_x, min_value_z
    REAL(DPR) :: max_value_x, max_value_z

    CALL print_general_msg('Update points boundary marker information')
    max_value_x = MAXVAL( nodeobj%x )
    max_value_z = MAXVAL( nodeobj%z )
    min_value_x = MinVAL( nodeobj%x )
    min_value_z = MinVAL( nodeobj%z )

    PRINT*, '------DEBUG---- Computational Domain infor ....'
    PRINT*, '------DEBUG---- min_value_x = ', REAL( MINVAL( nodeobj%x ) )
    PRINT*, '------DEBUG---- max_value_x = ', REAL( MAXVAL( nodeobj%x ) )
    PRINT*, '------DEBUG---- min_value_z = ', REAL( MINVAL( nodeobj%z ) )
    PRINT*, '------DEBUG---- max_value_z = ', REAL( MAXVAL( nodeobj%z ) )

    npoint = SIZE(nodeobj)
    DO i = 1, npoint
       ! check all other points with marker beyond 0 and 1 values
       IF( nodeobj(i)%bud /= 1 .AND. nodeobj(i)%bud /= 0 ) THEN
          IF( ABS( nodeobj(i)%x - max_value_x ) < tinyValue  .OR. &
               ABS( nodeobj(i)%x - min_value_x ) < tinyValue  .OR. &
               ABS( nodeobj(i)%z - max_value_z ) < tinyValue  .OR. &
               ABS( nodeobj(i)%z - min_value_z ) < tinyValue )  THEN

             nodeobj(i)%bud = 1
          end IF
       END IF
    end DO
  end SUBROUTINE check_node_boundary

  SUBROUTINE update_node_regional_marker(nregion, regionMark)

    !--set up regional marker for each node
    ! For a point shared by many cells:
    ! diff region's marker: 10   20  30  ... etc
    ! each element's reg:   1    0   0
    !                       0    1   0
    !                       0    1   0
    !                       0    0   1
    !                  -----------------
    !     in total          1    2   1
    ! which means: the point is shared by four elements from 3 diff regions
    !    with 1 ele from '10'-region, 
    !    with 2 eles from '20'-region, 
    !    with 1 ele from '30'-region, 
    !--Note: need mesh information: ele2node
    !        need the node entity
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nregion, regionMark(:)
    INTEGER :: npoint, nele
    INTEGER,ALLOCATABLE   :: nodeRegCounter(:,:)
    INTEGER    :: nnon
    INTEGER    :: i, j, k

    npoint = SIZE(nodes); nele = SIZE(ele2node, 1)
    IF( ALLOCATED(nodeRegCounter) .EQV. .TRUE.) DEALLOCATE(nodeRegCounter)
    ALLOCATE( nodeRegCounter(npoint, nregion) )
    nodeRegCounter = 0

    ele1: DO i = 1, nele
       DO j = 1, nregion
          ! regional marker of the element equals to which value in array of all markers?
          IF ( regionMark(j) == ele2node(i, 4) )  THEN
             !! 'j' is the right region index
             node3: DO k = 1, 3
                ! counter of the corresponding node pluses one
                nodeRegCounter( ele2node(i, k), j ) = nodeRegCounter( ele2node(i, k), j ) + 1
             end DO node3
          ENDIF
       ENDDO
    end DO ele1

    DO i = 1, npoint
       ! get the number of non-zero markers of all regions for each node
       nnon = COUNT( nodeRegCounter(i, :) /= 0)

       IF( nnon <= 0 )  THEN
          CALL print_debug_msg('At Point', i)
          CALL print_error_msg('Regional marker searching for this point failed')

       ELSEIF( nnon == 1 ) THEN
          ! when the node is shared by only one region,
          ! which region's marker for this node is nonzero?
          ! reg stores all the regional markers of those regions sharing the node
          ALLOCATE( nodes(i)%reg( nnon ) )
          nr:DO j = 1, nregion
             IF( nodeRegCounter(i, j) /= 0 )  THEN
                ! set up regional marker for this node
                nodes(i)%reg(1) =regionMark(j)
                EXIT nr
             ENDIF
          ENDDO nr
       ELSEIF( nnon >= 2 ) THEN
          ! when the node is shared by more than one regions
          ALLOCATE( nodes(i)%reg( nnon ) )
          nodes(i)%reg = 0
          k = 0
          DO j = 1, nregion
             IF( nodeRegCounter(i, j) /= 0 )  THEN
                ! set up reginal marker for the node
                ! -- store all the markers of those regions that share the node
                ! -- the markers are set by the user when meshing the domain, and have
                !    no repeated values (each marker is unique)
                k = k + 1
                nodes(i)%reg( k ) = regionMark(j)
             end IF
          end do
       ENDIF
    end DO
    DEALLOCATE( nodeRegCounter )
    RETURN
  end SUBROUTINE update_node_regional_marker



  ! ---------------------------------------  
  SUBROUTINE read_element_file(filepath, filename, attriMarker, eleobj)
    IMPLICIT NONE
    INTEGER :: fid, nele, i, inull
    CHARACTER(LEN=*),INTENT(IN) :: filepath, filename
    INTEGER,INTENT(IN) :: attriMarker(:)  ! elemental attribute markers
    INTEGER,INTENT(INOUT),ALLOCATABLE :: eleobj(:,:)

    fid = 41
    ! read .ele file
    OPEN(UNIT=fid, FILE= TRIM(ADJUSTL(filepath))//TRIM(ADJUSTL(filename)), STATUS = 'OLD')
    READ(fid, *) nele
    IF(ALLOCATED(eleobj) .EQV. .TRUE.)  DEALLOCATE( eleobj )
    ALLOCATE( eleobj(nele, 4) )
    DO i = 1, nele
       !--node1, node2, node3, ele-regional marker
       READ(fid, *) inull, eleobj(i,1), eleobj(i,2), eleobj(i,3), eleobj(i,4)

       ! If elemental attribute marker from ele file is outside the range of values that
       ! are defined when doing the meshing, then check such situations.
       IF( ANY( attriMarker == eleobj(i,4) ) .EQV. .FALSE. )  THEN
          PRINT*, ''; PRINT*, 'Error: the elemental region mark value: ', eleobj(i,4),' for the ', i, &
               '-th element/record from the file has no match in the defined mark values !'
          PRINT*, 'Program terminated by <STOP> !'
          STOP
       end IF
    end DO
    CLOSE(fid)
  end SUBROUTINE read_element_file
  ! ---------------------------------------
!!$  SUBROUTINE read_edge_file(filepath, filename, edgeobj)
!!$    IMPLICIT NONE
!!$    INTEGER :: fid, nedge, i, inull
!!$    CHARACTER(LEN=*),INTENT(IN) :: filepath, filename
!!$    INTEGER,INTENT(INOUT),ALLOCATABLE :: edgeobj(:,:)
!!$
!!$    fid = 41
!!$    ! read .edge file
!!$    OPEN(UNIT=fid, FILE= TRIM(ADJUSTL(filepath))//TRIM(ADJUSTL(filename)), STATUS = 'OLD')
!!$    READ(fid, *) nedge
!!$    ! read .edge file (as output file with Triangle (v1.6, 2005) flag '-e' used)
!!$    IF(ALLOCATED(edgeobj) .EQV. .TRUE.)  DEALLOCATE( edgeobj );  ALLOCATE( edgeobj(nedge, 3) )
!!$    DO i = 1, nedge
!!$       !--index of edge, node1, node2, bound marker
!!$       READ(fid, *) inull, edgeobj(i,1), edgeobj(i,2), edgeobj(i,3)
!!$    end DO
!!$    CLOSE(fid)
!!$  end SUBROUTINE read_edge_file

!!$  SUBROUTINE check_edge_boundary(nodeobj, edgeobj)
!!$    ! make sure the boundary edge markers are OK.
!!$    USE derived_data_module, ONLY: geoph_EM_node_type
!!$    IMPLICIT NONE
!!$    INTEGER :: i, nedge
!!$    TYPE(geoph_EM_node_type),INTENT(IN) :: nodeobj(:)
!!$    INTEGER,INTENT(INOUT),ALLOCATABLE :: edgeobj(:,:)
!!$
!!$    nedge = SIZE(edgeobj, 1)
!!$    DO i = 1, nedge
!!$       ! The edge boundary marker may not be reliable for it being used to detect actual boundary edges; instead,
!!$       !  the status of a boundary edge is calculated here: if both ends of an edge are boundary nodes, then the
!!$       !  edge is a boundary edge, otherwise it is not on the boundaries.
!!$       ! Best to be called/executed after reading .node file
!!$       IF( nodeobj(edgeobj(i,1))%bud == 1 .AND. nodeobj(edgeobj(i,2))%bud == 1 ) THEN
!!$          edgeobj(i,3) = 1  ! boundary edges
!!$       ELSE
!!$          edgeobj(i,3) = 0   ! non-boundary edges
!!$       ENDIF
!!$    end DO
!!$
!!$    RETURN
!!$  end SUBROUTINE check_edge_boundary

  ! ---------------------------------------
  SUBROUTINE read_neigh_element_file(filepath, filename, neiobj)
    ! neighbouring triangle facets for a particular triangle facet.
    IMPLICIT NONE
    INTEGER :: fid, i, inull, nele
    CHARACTER(LEN=*),INTENT(IN) :: filepath, filename
    INTEGER,INTENT(INOUT),ALLOCATABLE :: neiobj(:,:)

    fid = 41
    ! read .neigh file
    OPEN(UNIT=fid, FILE= TRIM(ADJUSTL(filepath))//TRIM(ADJUSTL(filename)), STATUS = 'OLD')
    READ(fid, *) nele
    IF(ALLOCATED(neiobj) .EQV. .TRUE.)  DEALLOCATE( neiobj );   ALLOCATE( neiobj(nele, 3) )
    DO i = 1, nele
       !--ele1, ele2, ele3 (-1:no neighbour for this side)
       READ(fid, *) inull, neiobj(i,1), neiobj(i,2), neiobj(i,3)
    end DO
    CLOSE(fid)
  end SUBROUTINE read_neigh_element_file

  ! ----------------------------------------------------------------------------------
  SUBROUTINE get_sub_nodes(mark, nodes, subnodes, node_g2l_map, nodemap)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: mark
    TYPE(node_cls),INTENT(INOUT),ALLOCATABLE :: nodes(:)
    TYPE(node_cls),INTENT(INOUT),ALLOCATABLE :: subnodes(:)
    INTEGER,INTENT(INOUT),ALLOCATABLE :: nodemap(:)  ! map from new sub-node ordering to the old ordering
    INTEGER,INTENT(INOUT),ALLOCATABLE :: node_g2l_map(:)
    INTEGER :: np, nreg, m, k

    np = SIZE(nodes)
    ALLOCATE( node_g2l_map(np))
    m = 0
    DO k = 1, np
       IF( ANY(nodes(k)%reg == mark) ) m = m + 1
    END DO
    ALLOCATE( subnodes(m))
    ALLOCATE( nodemap(m))

    m = 0
    node_g2l_map = 0
    nodemap = 0
    DO k = 1, np
       IF( ANY(nodes(k)%reg == mark) ) THEN
          m = m + 1
          subnodes(m)%x = nodes(k)%x
          subnodes(m)%z = nodes(k)%z
          subnodes(m)%bud = nodes(k)%bud
          nreg = SIZE(nodes(k)%reg)
          IF(.NOT. ALLOCATED(subnodes(m)%reg) ) ALLOCATE(subnodes(m)%reg(nreg))
          subnodes(m)%reg = nodes(k)%reg

          ! update sub-node boundary info
          IF(nreg >= 2) subnodes(m)%bud = 1

          nodemap(m) = k  ! original global index of this node

          node_g2l_map(k) = m
       END IF
    END DO
    RETURN
  end SUBROUTINE get_sub_nodes

  ! ----------------------------------------------------------------------------------
  SUBROUTINE get_sub_elements(mark, ele2node, neighele, node_g2l_map, new_ele2node, new_neighele)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: mark, ele2node(:,:), neighele(:,:), node_g2l_map(:)
    INTEGER,INTENT(INOUT),ALLOCATABLE :: new_ele2node(:,:), new_neighele(:,:)
    INTEGER,ALLOCATABLE :: ele_g2l_map(:)  ! map from  the old ordering to the new sub-ele ordering
    INTEGER :: ne, m, k, j, wp

    ne = SIZE(ele2node, 1)
    ALLOCATE( ele_g2l_map(ne) )
    m = 0
    DO k = 1, ne
       IF( ele2node(k, 4) == mark ) m = m + 1
    END DO
    ALLOCATE( new_ele2node(m, 4), new_neighele(m, 3))

    m = 0
    ele_g2l_map = 0  ! global index of element to local index map

    DO k = 1, ne
       IF( ele2node(k, 4) == mark ) THEN
          m = m + 1
          new_ele2node(m, :) = ele2node(k, :) ! for now, just copy the global 3 nodal indexes as well

          new_neighele(m,:) = neighele(k,:) ! for now, just copy the global neighbouring elements
          !elemap(m) = k  ! original global index of this element

          ele_g2l_map(k) = m
       END IF
    END DO

    ! update nodal indexes in the new ele2node
    ne = SIZE(new_ele2node, 1)
    DO k = 1, ne
       DO j = 1, 3
          wp = new_ele2node(k, j)
          new_ele2node(k, j) = node_g2l_map(wp)
          IF( node_g2l_map(wp) ==0 ) THEN
             PRINT*, 'The global node is ', wp
             CALL print_error_msg('Global-to-local node map is 0 !')
          END IF
       END DO
    end DO

    ! update elemental indexes in the new neighele
    DO k = 1, ne
       in: DO j = 1, 3   ! maximum 3 neighbouring elements
          wp = new_neighele(k, j)  ! original elemental index
          IF( wp == -1 ) CYCLE in
          IF( ele2node(wp, 4) == mark ) THEN
             ! this neighbour ele is in the considered region
             new_neighele(k, j) = ele_g2l_map(wp)
          ELSE
             ! the neighbour ele is outside the region
             new_neighele(k, j) = -1
          END IF
       end DO in
    END DO
    RETURN
  end SUBROUTINE get_sub_elements
    ! ---------------------------------------

  SUBROUTINE write_mesh_files(n, mesh_path, basename, new_nodes, new_ele2node, new_neighele, nodemap, node_g2l_map)
    IMPLICIT NONE

    INTEGER :: fid, np, ne, k
    CHARACTER(LEN=*),INTENT(IN) :: mesh_path, basename
    TYPE(node_cls),INTENT(IN) :: new_nodes(:)
    INTEGER,INTENT(IN) :: n, new_ele2node(:,:), new_neighele(:,:), nodemap(:), node_g2l_map(:)
    CHARACTER(LEN=400) :: nodefile, elefile, neifile, mapfile, imapfile
    CHARACTER(LEN=2) :: sn

    WRITE(sn,'(I2.2)') n
    nodefile = TRIM(ADJUSTL(mesh_path))//TRIM(ADJUSTL(basename))//".sub"//sn//".node"
    elefile = TRIM(ADJUSTL(mesh_path))//TRIM(ADJUSTL(basename))//".sub"//sn//".ele"
    neifile = TRIM(ADJUSTL(mesh_path))//TRIM(ADJUSTL(basename))//".sub"//sn//".neigh"
    mapfile = TRIM(ADJUSTL(mesh_path))//TRIM(ADJUSTL(basename))//".sub"//sn//".map"
    imapfile = TRIM(ADJUSTL(mesh_path))//TRIM(ADJUSTL(basename))//".sub"//sn//".invmap"
    PRINT*, "Writing mesh files:"
    PRINT*, TRIM(ADJUSTL(nodefile))
    PRINT*, TRIM(ADJUSTL(elefile))
    PRINT*, TRIM(ADJUSTL(neifile))
    fid = 40
    ! write .node file
    OPEN(UNIT=fid, FILE= nodefile, STATUS = 'unknown')
    np = SIZE(new_nodes)
    WRITE(fid,*) np, 2, 0, 1
    DO k = 1, np
       ! In 2-D mesh, y is undefined.
       WRITE(fid,*)  k, new_nodes(k)%x, new_nodes(k)%z, new_nodes(k)%bud
    END DO
    CLOSE(fid)

    ! .ele file
    OPEN(UNIT=fid, FILE= elefile, STATUS = 'unknown')
    ne = SIZE(new_ele2node, 1)
    WRITE(fid,*) ne, 3, 1
    DO k = 1, ne
       WRITE(fid,*) k, new_ele2node(k,:)
    END DO
    CLOSE(fid)

    ! .neigh file
    OPEN(UNIT=fid, FILE= neifile, STATUS = 'unknown')
    ne = SIZE(new_neighele, 1)
    WRITE(fid,*) ne, 3
    DO k = 1, ne
       WRITE(fid,*) k, new_neighele(k,:)
    END DO
    CLOSE(fid)

    ! .map file
    OPEN(UNIT=fid, FILE= mapfile, STATUS = 'unknown')
    np = SIZE(nodemap)
    WRITE(fid,*) np
    DO k = 1, np
       WRITE(fid,*) k, nodemap(k)
    END DO
    CLOSE(fid)

        ! .invmap file
    OPEN(UNIT=fid, FILE= imapfile, STATUS = 'unknown')
    np = SIZE(node_g2l_map)
    WRITE(fid,*) np
    DO k = 1, np
       WRITE(fid,*) k, node_g2l_map(k)
    END DO
    CLOSE(fid)

  end SUBROUTINE write_mesh_files

  !! ------------------------------------------------------------------------------------------
end PROGRAM mesh2d_partition





