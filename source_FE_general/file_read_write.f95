
MODULE file_rw
  ! Author: Jianbo Long
  USE float_precision, ONLY: DPR, CDPR
  USE derived_data_module, ONLY: nodes
  USE modelling_parameter, ONLY: modelling, EM_generic
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: write_to_files_Ron_MT2D, read_sol_Ron_MT2D, &
       write_numerical_result_to_files

  INTERFACE write_numerical_result_to_files
     MODULE PROCEDURE write_complex_numerical_result_to_files, write_real_numerical_result_to_files
  END INTERFACE write_numerical_result_to_files
CONTAINS
  SUBROUTINE write_to_files_Ron_MT2D(solEy, solHy)
    ! for domain decomposition tests
    USE modelling_parameter, ONLY: mesh_par
    IMPLICIT NONE
    INTEGER :: nobs
    REAL(DPR),ALLOCATABLE :: xobs(:), zobs(:)
    COMPLEX(CDPR), INTENT(IN) :: solEy(:), solHy(:)
    INTEGER :: fid, k, j, regionMark
    CHARACTER(LEN=2) :: nsub
    CHARACTER(:),ALLOCATABLE :: str, sufix

    fid = 48
    nobs = SIZE(nodes)
    ALLOCATE(xobs(nobs), zobs(nobs))
    xobs = nodes%x
    zobs = nodes%z
    str = TRIM(ADJUSTL(mesh_par%basefilename))
    k = LEN( str )
    IF( INDEX(str, "sub") > 0 .OR. INDEX(str, "SUB") > 0 ) THEN
       ! if 'str' contains "sub"
       nsub = str(k-1: k)
       sufix = "sol_"//nsub
    ELSE
       sufix = "sol"
    END IF

    ! TE
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//'MT2D_TE_'//sufix//'.txt')
    WRITE(fid, *) 'Data NO. & Freq(Hz)',  ' & ', 'position-m ( x, z )',  ' || ', 'sol-Real',  ' & ','sol-Img', &
         ' || ', 'RegionMark',  ' & ','BoundaryMark'
    WRITE(fid, *) nobs
    DO k = 1, nobs
       IF( SIZE(nodes(k)%reg) <= 1 ) THEN
          regionMark = nodes(k)%reg(1)
       ELSE
          regionMark = 1
          DO j = 1, SIZE(nodes(k)%reg)
             regionMark = regionMark * nodes(k)%reg(j)
          END DO
       END IF
       WRITE(fid, '(I6, G14.5, 2F14.3, 2G28.16, 2I9)') k, EM_generic%freq, xobs(k), zobs(k), &
            REAL(solEy(k), kind=8), AIMAG(solEy(k) ), regionMark, nodes(k)%bud
    end DO
    CLOSE(UNIT=fid)

    ! TM 
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//'MT2D_TM_'//sufix//'.txt')
    WRITE(fid, *) 'Data NO. & Freq(Hz)',  ' & ', 'position-m ( x, z )',  ' || ', 'sol-Real',  ' & ','sol-Img', &
         ' || ', 'RegionMark',  ' & ','BoundaryMark'
    WRITE(fid, *) nobs
    DO k = 1, nobs
       IF( SIZE(nodes(k)%reg) <= 1 ) THEN
          regionMark = nodes(k)%reg(1)
       ELSE
          regionMark = 1
          DO j = 1, SIZE(nodes(k)%reg)
             regionMark = regionMark * nodes(k)%reg(j)
          END DO
       END IF
       WRITE(fid, '(I6, G14.5, 2F14.3, 2G28.16, 2I9)') k, EM_generic%freq, xobs(k), zobs(k), &
            REAL(solHy(k), kind=8), AIMAG(solHy(k) ), regionMark, nodes(k)%bud
    end DO
    CLOSE(UNIT=fid)

    RETURN
  END SUBROUTINE write_to_files_Ron_MT2D

  SUBROUTINE read_sol_Ron_MT2D(solEy, solHy, x, z, nodemap)
    ! for Ron's domain decomposition tests -- read global solutions for getting sub-domain solution
    USE modelling_parameter, ONLY: mesh_par
    IMPLICIT NONE
    INTEGER,INTENT(INOUT),ALLOCATABLE :: nodemap(:)
    REAL(DPR),INTENT(INOUT),ALLOCATABLE :: x(:), z(:)
    COMPLEX(CDPR),INTENT(INOUT),ALLOCATABLE :: solEy(:), solHy(:)
    INTEGER :: fid, k, ndata, null, nu2, nu3
    REAL(DPR) :: freq, real_p, imag_p

    fid = 48

    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//'MT2D_TE_sol.txt', STATUS = 'OLD')
    READ(fid, *)    ! skip header
    READ(fid, *) ndata
    ALLOCATE(x(ndata), z(ndata), solEy(ndata), solHy(ndata))
    DO k = 1, ndata
       READ(fid, '(I6, G14.5, 2F14.3, 2G28.16, 2I9)') null, freq, x(k), z(k), &
            real_p, imag_p, nu2, nu3
       solEy(k) = COMPLEX(real_p, imag_p)
    end DO
    CLOSE(UNIT=fid)

    ! for debugging
    IF( ABS(EM_generic%freq - freq)/ EM_generic%freq  >= 0.1) THEN
       PRINT*, 'Freqs not matching !!'
       PRINT*, "EM_generic%freq = ", EM_generic%freq
       PRINT*, "freq = ", freq
       STOP
    END IF

    ! TM 
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//'MT2D_TM_sol.txt', STATUS = 'OLD')
    READ(fid, *)    ! skip header
    READ(fid, *) ndata
    DO k = 1, ndata
       READ(fid, '(I6, G14.5, 2F14.3, 2G28.16, 2I9)') null, freq, x(k), z(k), &
            real_p, imag_p, nu2, nu3
       solHy(k) = COMPLEX(real_p, imag_p)
    end DO
    CLOSE(UNIT=fid)

    ! .map file
    OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(mesh_par%meshfilepath))//TRIM(ADJUSTL(mesh_par%basefilename))//".map", STATUS = 'OLD')
    !READ(fid, *)    ! skip header
    READ(fid, *) ndata
    ALLOCATE( nodemap(ndata) )
    DO k = 1, ndata
       READ(fid, *) null, nodemap(k)
    end DO
    CLOSE(UNIT=fid)

    RETURN
  END SUBROUTINE read_sol_Ron_MT2D

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------
  SUBROUTINE read_surface_mapping_RHS_values()
    IMPLICIT NONE
    RETURN
  end SUBROUTINE read_surface_mapping_RHS_values

  ! ----------------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------------
  SUBROUTINE write_complex_numerical_result_to_files(nobs, xobs, yobs, zobs, func, polar)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nobs
    REAL(DPR),INTENT(IN) :: xobs(nobs), yobs(nobs), zobs(nobs)
    COMPLEX(CDPR), INTENT(IN) :: func(nobs)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: polar
    INTEGER :: fid, k

    fid = 41
    IF( PRESENT(polar) ) THEN
       OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//&
            'FE_'//TRIM(ADJUSTL(polar))//'_result.txt')
    ELSE
       OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//'FE_result.txt')
    END IF

    WRITE(fid, *) 'Data NO.',  ' & ', 'Pos ( x, y, z )',  ' || ', 'F--Real',  ' & ','F--Img'
    DO k = 1, nobs
       WRITE(fid, '(I6, 3F14.3, 2G18.7)') k, xobs(k), yobs(k), zobs(k), REAL(func(k), kind=8), AIMAG(func(k) )
    end DO
    CLOSE(UNIT=fid)
    RETURN
  END SUBROUTINE write_complex_numerical_result_to_files
  ! ----------------------------------------------------------------------------------------
  SUBROUTINE write_real_numerical_result_to_files(nobs, xobs, yobs, zobs, func, polar)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nobs
    REAL(DPR),INTENT(IN) :: xobs(nobs), yobs(nobs), zobs(nobs)
    REAL(DPR), INTENT(IN) :: func(nobs)
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: polar
    INTEGER :: fid, k

    fid = 41
    IF( PRESENT(polar) ) THEN
       OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//&
            'FE_'//TRIM(ADJUSTL(polar))//'_result.txt')
    ELSE
       OPEN(UNIT=fid, FILE=TRIM(ADJUSTL(modelling%outputfilepath))//'FE_result.txt')
    END IF

    WRITE(fid, *) 'Data NO.',  ' & ', 'Pos ( x, y, z )',  ' || ', 'Func'
    DO k = 1, nobs
       WRITE(fid, '(I6, 3F14.3, G18.7)') k, xobs(k), yobs(k), zobs(k), func(k)
    end DO
    CLOSE(UNIT=fid)
    RETURN
  END SUBROUTINE write_real_numerical_result_to_files

end MODULE file_rw
