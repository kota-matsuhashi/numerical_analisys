program report_rocket3
    implicit none

    integer, parameter :: imax=21, jmax=21, kmax=21
    integer, parameter :: nmax=30000
    integer, parameter :: case=1     !計算条件(Case 1 or 2 or 3)
    integer :: i, j, k, n
    real, parameter :: dt=1.0d0, dx=1.0d-3
    real, allocatable, dimension(:,:,:) :: T, T_new
    real, allocatable, dimension(:,:,:) :: lambda, rho, c
    real :: T0, T1, T2, T3, T4, T5, T6
    real :: T_sum, T_ave
    real, parameter :: T_a=2.0d1
    real, parameter :: beta_l=1.0d0
    real, parameter :: rho_l=1.0d3, rho_s=9.2d2
    real, parameter :: c_l=4.2d3
    real, parameter :: lambda_l=5.7d-1, lambda_s=2.2d0
    real :: lambda_m

    allocate(T(imax,jmax,kmax))
    allocate(T_new(imax,jmax,kmax))
    allocate(lambda(imax,jmax,kmax))
    allocate(rho(imax,jmax,kmax))
    allocate(c(imax,jmax,kmax))

    ! 初期条件設定
    call initialize(T, lambda, rho, c)

    ! 時間ステップループ
    do n = 1, nmax
        call update_temperature(T, T_new, lambda, rho, c, dt, dx)
        T = T_new
    end do

    ! 結果の出力
    call output_results(T)

    deallocate(T)
    deallocate(T_new)
    deallocate(lambda)
    deallocate(rho)
    deallocate(c)

end program report_rocket3

subroutine initialize(T, lambda, rho, c)
    implicit none
    integer, parameter :: imax=21, jmax=21, kmax=21
    real, dimension(imax,jmax,kmax) :: T, lambda, rho, c
    integer :: i, j, k

    ! 初期温度設定
    T = 20.0

    ! 物性値設定
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                if (k <= kmax/2) then
                    lambda(i,j,k) = 5.7d-1
                    rho(i,j,k) = 1.0d3
                    c(i,j,k) = 4.2d3
                else
                    lambda(i,j,k) = 2.2d0
                    rho(i,j,k) = 9.2d2
                    c(i,j,k) = 2.0d3
                end if
            end do
        end do
    end do
end subroutine initialize

subroutine update_temperature(T, T_new, lambda, rho, c, dt, dx)
    implicit none
    integer, parameter :: imax=21, jmax=21, kmax=21
    real, dimension(imax,jmax,kmax) :: T, T_new, lambda, rho, c
    real :: dt, dx
    integer :: i, j, k
    real :: alpha

    do k = 2, kmax-1
        do j = 2, jmax-1
            do i = 2, imax-1
                alpha = lambda(i,j,k) / (rho(i,j,k) * c(i,j,k))
                T_new(i,j,k) = T(i,j,k) + alpha * dt / dx**2 * &
                               (T(i+1,j,k) - 2.0*T(i,j,k) + T(i-1,j,k) + &
                                T(i,j+1,k) - 2.0*T(i,j,k) + T(i,j-1,k) + &
                                T(i,j,k+1) - 2.0*T(i,j,k) + T(i,j,k-1))
            end do
        end do
    end do
end subroutine update_temperature

subroutine output_results(T)
    implicit none
    integer, parameter :: imax=21, jmax=21, kmax=21
    real, dimension(imax,jmax,kmax) :: T
    integer :: i, j, k

    open(unit=10, file='temperature_results.dat', status='unknown')
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                write(10,*) i, j, k, T(i,j,k)
            end do
        end do
    end do
    close(10)
end subroutine output_results