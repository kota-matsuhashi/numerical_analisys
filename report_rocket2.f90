program report_rocket2
    implicit none

    integer,parameter :: case=1     !計算条件(Case 1 or 2 or 3)

    integer n
    integer, parameter :: nmax=30000
    integer i, j, k
    integer, parameter :: imax=21, jmax=21, kmax=21
    integer count_w
    real,parameter :: dt=1.0d0, dx=1.0d-3
    real,allocatable,dimension(:,:,:) :: T, T_re
    real,allocatable,dimension(:,:,:) :: beta, beta_re
    real T0, T1, T2, T3, T4, T5, T6
    real bt
    real T_sum, T_ave
    real,parameter :: T_a=2.0d1
    real,parameter :: beta_l=1.0d0
    real,parameter :: rho_l=1.0d3, rho_s=9.2d2
    real,parameter :: c_l=4.2d3
    real,parameter :: lambda_l=5.7d-1, lambda_s=2.2d0
    real lambda_m
    real alpha
    real mu
    real,parameter :: q_w=3.3d5
    real phi, omega
    real e1, t_1(2), e2, t_2(2)

    e1 = etime(t_1)

    allocate(T(imax, jmax, kmax))
    allocate(beta(imax, jmax, kmax))
    allocate(T_re(imax, jmax, kmax))
    allocate(beta_re(imax, jmax, kmax))

    !初期条件
    call init(case, T(:,:,:), beta(:,:,:), imax, jmax, kmax)

    !水面の境界条件
    do j = 1, jmax
        do i = 1, imax
            T(i, j, kmax) = T_a
        enddo
    enddo

    !初期値の可視化ファイル出力．このときに格子データも出力しておく(mode=1)
    call graphic(case, 0, T(:,:,:), beta(:,:,:), imax, jmax, kmax, dx)

    alpha = lambda_l / (rho_l * c_l)
    mu = alpha * dt / dx **2

    if(case .eq. 1) then
        phi = (2.0d1 / 2.1d1) ** 3
    else if(case .eq. 2) then
        phi = (1.0d1 / 1.1d1) ** 3
    endif

    open(10, file='T_average.dat', status='replace')
    write(10, 50) 'second', 'T_ave'
    close(10)

    do n = 1, nmax

        do k = 1, kmax-1
            do j = 1, jmax
                do i = 1, imax
                    if(beta(i, j, k) .eq. beta_l) then  !水の場合
                        T0 = T(i, j, k)

                        if(i .eq. 1) then
                            T1 = T(i+1, j, k)
                            T2 = T(i+1, j, k)
                        else if(i .eq. imax) then
                            T1 = T(i-1, j, k)
                            T2 = T(i-1, j, k)
                        else
                            T1 = T(i+1, j, k)
                            T2 = T(i-1, j, k)
                        endif

                        if(j .eq. 1) then
                            T3 = T(i, j+1, k)
                            T4 = T(i, j+1, k)
                        else if(j .eq. jmax) then
                            T3 = T(i, j-1, k)
                            T4 = T(i, j-1, k)
                        else
                            T3 = T(i, j+1, k)
                            T4 = T(i, j-1, k)
                        endif

                        if(k .eq. 1) then
                            T5 = T(i, j, k+1)
                            T6 = T(i, j, k+1)
                        else
                            T5 = T(i, j, k+1)
                            T6 = T(i, j, k-1)
                        endif

                        !熱伝導方程式
                        T_re(i, j, k) = T0 + mu*(T1+T2+T3+T4+T5+T6-6.0d0*T0)
                    else                    !氷の場合
                        bt = beta(i, j, k)
                        lambda_m = ((1+bt)*lambda_l + (1-bt)*lambda_s) / 2.0d0
                        omega = lambda_m * dt / (phi * rho_s * q_w * dx **2)
                        T0 = T(i, j, k)
                        T1 = T(i+1, j, k)
                        T2 = T(i-1, j, k)
                        T3 = T(i, j+1, k)
                        T4 = T(i, j-1, k)
                        T5 = T(i, j, k+1)
                        T6 = T(i, j, k-1)

                        !β(水の質量分率)を更新
                        beta_re(i, j, k) = bt + omega*(T1+T2+T3+T4+T5+T6-6.0d0*T0)

                        !βが1を超えた点はその後水として扱う
                        if(beta_re(i, j, k) .ge. beta_l) then
                            beta_re(i, j, k) = beta_l
                        endif
                    endif
                enddo
            enddo
        enddo

        T_sum = 0.0d0
        count_w = 0
        do k = 1, kmax-1
            do j = 1, jmax
                do i = 1, imax
                    if(beta(i, j, k) .eq. beta_l) then
                        T(i, j, k) = T_re(i, j, k)
                    else
                        beta(i, j, k) = beta_re(i, j, k)
                    endif

                    if(beta(i, j, k) .eq. beta_l) then
                        T_sum = T_sum + T(i, j, k)
                        count_w = count_w + 1
                    endif
                enddo
            enddo
        enddo

        !水の平均温度を出力(50ステップに1回)
        if(mod(n,50) .eq. 0) then
            T_ave = T_sum / count_w
            open(10, file='T_average.dat', access='append')
            write(10, 60) n, T_ave
            close(10)
        endif

        !可視化ファイルの出力(500ステップに1回)．格子データは出力しない(mode=2)
        if(mod(n, 500) .eq. 0 .and. n .le. 3000) then
            call graphic(2, n, T(:,:,:), beta(:,:,:), imax, jmax, kmax, dx)
        endif

    enddo

    e2 = etime(t_2)
    write(*, *) 'elapsed time:', e2-e1
    write(*, *) 'user:', t_2(1)-t_1(1), 'system:', t_2(2)-t_1(2)

    50 format(a6, 1x, a5)
    60 format(i6, 1x, f6.3)

end program report_rocket2

subroutine init(case, T, beta, imax, jmax, kmax)
    integer case
    integer i, j, k
    integer imax, jmax, kmax
    integer ii, jj, kk
    real,dimension(imax,jmax,kmax) :: T, beta
    real,parameter :: T_l=5.0d0, T_s=0.0d0
    real,parameter :: beta_l=1.0d0, beta_s=0.0d0
    
    do k = 1, kmax
        do j = 1, jmax
            do i = 1, imax
                T(i, j, k) = T_l
                beta(i, j, k) = beta_l
            enddo
        enddo
    enddo

    if(case .eq. 1) then
        ii = int(imax/3)
        jj = int(jmax/3)
        kk = int(kmax/3)
        do k = kk+1, 2*kk+1
            do j = jj+1, 2*jj+1
                do i = ii+1, 2*ii+1
                    T(i, j, k) = T_s
                    beta(i, j, k) = beta_s
                enddo
            enddo
        enddo
    else if(case .eq. 2) then
        ii = int(imax/6)
        jj = int(jmax/6)
        kk = int(kmax/6)
        do k = kk+1, 2*kk+1
            do j = jj+1, 2*jj+1
                do i = ii+1, 2*ii+1
                    T(i, j, k) = T_s
                    beta(i, j, k) = beta_s
                    T(imax-i+1, j, k) = T_s
                    beta(imax-i+1, j, k) = beta_s
                    T(i, jmax-j+1, k) = T_s
                    beta(i, jmax-j+1, k) = beta_s
                    T(imax-i+1, jmax-j+1, k) = T_s
                    beta(imax-i+1, jmax-j+1, k) = beta_s
                    T(i, j, kmax-k+1) = T_s
                    beta(i, j, kmax-k+1) = beta_s
                    T(imax-i+1, j, kmax-k+1) = T_s
                    beta(imax-i+1, j, kmax-k+1) = beta_s
                    T(i, jmax-j+1, kmax-k+1) = T_s
                    beta(i, jmax-j+1, kmax-k+1) = beta_s
                    T(imax-i+1, jmax-j+1, kmax-k+1) = T_s
                    beta(imax-i+1, jmax-j+1, kmax-k+1) = beta_s
                enddo
            enddo
        enddo
    endif

    return

end subroutine init

subroutine graphic(mode, step, T, beta, imax, jmax, kmax, dx)
    integer mode
    integer step
    integer i, j, k
    integer imax, jmax, kmax
    real dx
    real,dimension(imax,jmax,kmax) :: x, y, z
    real,dimension(imax,jmax,kmax) :: T, beta
    character*20 rsltname

    if(mode .eq. 1) then
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    x(i, j, k) = (i - 1) * dx
                    y(i, j, k) = (j - 1) * dx
                    z(i, j, k) = (k - 1) * dx
                enddo
            enddo
        enddo

        open(11, file='grid_motor.DAT', status='unknown', blank='null')
            write(11, 80) imax, jmax, kmax
            write(11, 90) real(x(1:imax,1:jmax,1:kmax))
            write(11, 90) real(y(1:imax,1:jmax,1:kmax))
            write(11, 90) real(z(1:imax,1:jmax,1:kmax))
        close(11)
    end if

    write(rsltname, 100) step

    open(12, file=rsltname, status='unknown', blank='null')
    write(12, 80) imax, jmax, kmax
    write(12, 80) 2
    write(12, 90) real(T(1:imax,1:jmax,1:kmax))
    write(12, 90) real(beta(1:imax,1:jmax,1:kmax))
    close(12)

    80 format(7i9)
    90 format(7f10.5)
    100 format('RSLT',i5.5,'.P3D')

    return
end subroutine graphic