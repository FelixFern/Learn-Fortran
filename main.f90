PROGRAM SWE
    IMPLICIT NONE

    INTEGER :: count_rate

    INTEGER, PARAMETER :: L = 25
    REAL, PARAMETER :: dx = 0.2
    INTEGER, PARAMETER :: T = 20000
    INTEGER, PARAMETER :: d0 = 10
    REAL, PARAMETER :: GRAVITY = 9.81
    INTEGER :: start_time, end_time

    INTEGER :: Nx, Nt
    REAL :: dt
    
    REAL :: elapsed_time
    INTEGER :: i, j
    REAL, DIMENSION(:), ALLOCATABLE :: L_arr, T_arr
    REAL, DIMENSION(:,:), ALLOCATABLE :: u, eta
    
    CALL SYSTEM_CLOCK(COUNT_RATE=count_rate)
    CALL SYSTEM_CLOCK(COUNT=start_time)

    Nx = INT(FLOOR(L / dx) + 1)
    ALLOCATE(L_arr(0:Nx))
    CALL linspace(0.0, REAL(L), L_arr)
    
    dt = dx / SQRT(GRAVITY * REAL(d0))
    Nt = INT(FLOOR(T / dt) + 1)
    ALLOCATE(T_arr(0:Nt))
    CALL linspace(0.0, REAL(T), T_arr)
    
    ALLOCATE(u(0:Nx, 0:Nt))
    ALLOCATE(eta(0:Nx, 0:Nt))
    
    u(:, 1) = 0.0
    DO i = 0, Nx
      eta(i, 1) = COS(3.14 / REAL(L) * L_arr(i))
    END DO
    
    DO i = 1, Nt - 1
      DO j = 0, Nx
        eta(j, i + 1) = eta(j, i) - dt / dx * REAL(d0) * (u(j + 1, i) - u(j, i))
      END DO
      
      DO j = 1, Nx
        u(j, i + 1) = u(j, i) - dt / dx * GRAVITY * (eta(j, i + 1) - eta(j - 1, i + 1))
      END DO
    END DO
    
    OPEN(UNIT=10, FILE='eta_data.txt', STATUS='replace', ACTION='write', FORM='formatted')
    
    DO i = 0, Nt
        WRITE(10, *) T_arr(i), eta(:, i)
    END DO
    
    CLOSE(UNIT=10)
    DEALLOCATE(L_arr, T_arr, u, eta)
    
    CALL SYSTEM_CLOCK(COUNT=end_time)

    elapsed_time = (end_time - start_time) / count_rate

    PRINT *, 'Computing time:', elapsed_time, ' seconds'
    
    PRINT *, Nx, Nt

    CONTAINS
        SUBROUTINE linspace(start, end, arr)
            REAL, INTENT(IN) :: start, end
            REAL, DIMENSION(:), INTENT(OUT) :: arr
            INTEGER :: i, n
            
            n = SIZE(arr) - 1
            DO i = 1, n
                arr(i) = start + REAL(i) / REAL(n) * (end - start)
            END DO

        END SUBROUTINE linspace

END PROGRAM SWE
  