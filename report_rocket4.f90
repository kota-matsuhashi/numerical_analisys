module initialize_params
    implicit none

    integer, parameter :: if=21, jf=21, kf=21
    integer :: ia, ib, ja, jb, ka, kb, na
    common /grid_params/ dx, dy, dz
    real :: dx, dy, dz, dx2, dy2, dz2, c0, omega, ck1, ck2
    common, real, dimension(if, jf, kf) :: p, u, v, w, T, al

    ! 定数初期化
    omega = 1.5
    ia, ib = 8, 14
    ja, jb = 8, 14
    ka, kb = 8, 14
    ck1, ck2 = 10.0, 1.0

end module initialize_params

program report_rocket4
    use initialize_params
    implicit none

    integer :: i, j, k, kp1, km1
    real :: dz1

    ! 初期化
    call initialize()

    ! 境界条件
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                T(i, j, k) = 50.0 ! 初期温度
                al(i, j, k) = 50.0 
            end do
        end do
    end do
    do k = 1, kf
        do j = 1, jf
            T(1, j, k) = 100.0 ! x=0での温度
            T(if, j, k) = 50.0 ! x=Lでの温度
        end do
    end do

    ! 計算ループ
    iter = 0
    tsum = 0.0
    iter = iter + 1
    tsum = tsum + dt

    do i = 2, if-1
        do j = 1, jf
            do k = 1, kf
                ! x方向
                im1 = i - 1
                if (i == 1) im1 = 2
                ip1 = i + 1
                if (i == if) im1 = if - 1                
                ckxl = ck1
                ckxr = ck1
                if ((j >= ja) .and. (j <= jb) .and. (k >= ka) .and. (k <= kb) .and. (i > ia) .and. (i <= ib)) ckxl = ck2
                if ((j >= ja) .and. (j <= jb) .and. (k >= ka) .and. (k <= kb) .and. (i >= ia) .and. (i < ib)) ckxr = ck2
                ! y方向
                km1 = k - 1
                if (k == 1) km1 = 2
                kp1 = k + 1
                if (k == kf) kp1 = kf - 1
                ckzl = ck1
                ckzr = ck1
                if ((j >= ja) .and. (j <= jb) .and. (k > ka) .and. (k <= kb) .and. (i >= ia) .and. (i <= ib)) ckzl = ck2
                if ((j >= ja) .and. (j <= jb) .and. (k >= ka) .and. (k < kb) .and. (i >= ia) .and. (i <= ib)) ckzr = ck2

                r = (ckxl * (T(ip1, j, k) - T(i, j, k)) - ckxr * (T(i, j, k) - T(im1, j, k))) * dlx * 0.5 &
                    + (ckyl * (T(i, jp1, k) - T(i, j, k)) - ckyr * (T(i, j, k) - T(i, jm1, k))) * dly * 0.5 &
                    + (ckzl * (T(i, j, kp1) - T(i, j, k)) - ckzr * (T(i, j, k) - T(i, j, km1))) * dlz * 0.5

                al(i, j, k) = r
            end do
        end do
    end do

    na = 0
    na = na + 1
    rmax = 0.0

    do ii = 2, if-1
        i = ii
        if(mod(na, 2) == 0) i = if - i + 1
        do jj = 1, jf
            j = jj
            if(mod(na, 2) == 0) j = jf - j + 1
            do kk = 1, kf
                k = kk
                if(mod(na, 2) == 0) k = kf - k + 1
                    
                dl = 1.0 + (dlx * (ckxr + ckxl) + dly * (ckyr + ckyl) + dlz * (ckzr + ckzl)) * 0.5
                r = (ckxl * T(im1, j, k) + ckxr * T(ip1, j, k)) * dlx / (dl * 2.0) &
                    + (ckyl * T(i, jm1, k) + ckyr * T(i, jp1, k)) * dly / (dl * 2.0) &
                    + (ckzl * T(i, j, km1) + ckzr * T(i, j, kp1)) * dlz / (dl * 2.0) &
                    + al(i, j, k) / dl - T(i, j, k)
                T(i, j, k) = T(i, j, k) + omega * r
                rmax = max(rmax, abs(r))
            end do
        end do
    end do

    if (na == 1) rmax1 = rmax
    rmax = rmax / rmax1
    write(6, 20) na, rmax
    write(11, 20) na, rmax
    format(5x, 'na = ',i3,5x,' rmax = ', f10.7)
    
    ! 計算データの後処理
    call graphic()

    stop
end program report_rocket4

subroutine initialize()
    use grid_params
    implicit none

    integer :: i, j, k
    ! 格子間隔
    dx = 0.1
    dy = 0.1
    dz = 0.1
    dx2 = dx * dx
    dy2 = dy * dy
    dz2 = dz * dz
    dlx = dt / dx2
    dly = dt / dy2
    dlz = dt / dz2

    ! 初期化
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                p(i, j, k) = 0.0
                u(i, j, k) = 0.0
                v(i, j, k) = 0.0
                w(i, j, k) = 0.0
                al(i, j, k) = 50.0
            end do
        end do
    end do

end subroutine initialize

subroutine graphic()
    use initialize_params
    implicit none

    integer :: i, j, k
    real, dimension(if) :: x
    real, dimension(jf) :: y
    real, dimension(kf) :: z

    ! dx, dy, dz の出力
    print *, 'dx = ', dx
    print *, 'dy = ', dy
    print *, 'dz = ', dz

    ! グリッドデータ出力
    x(1) = 0.0
    do i = 2, if
        x(i) = x(i-1) + dx
    end do
    y(1) = 0.0
    do j = 2, jf
        y(j) = y(j-1) + dy
    end do
    z(1) = 0.0
    do k = 2, kf
        z(k) = z(k-1) + dz
    end do

    ! グリッドデータ書き出し
    open(unit=11, file='report_grid.DAT', status='unknown', blank='null')
    write(11,*) if, jf, kf
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                write(11,*) x(i, j, k)
            end do
        end do
    end do
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                write(11,*) y(i, j, k)
            end do
        end do
    end do
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                write(11,*) z(i, j, k)
            end do
        end do
    end do
    close(11)

    ! 計算データ書き出し
    open(unit=12, file='report_data.DAT', status='unknown', blank='null')
    uin=1.0
    vis=0.0
    tin=0.0
    time=0.0
    write(12,*) if, jf, kf
    write(12,*) uin, vis, tin, time
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                write(12,*) p(i, j, k)!圧力
            end do
        end do
    end do
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                write(12,*) u(i, j, k)!x方向速度
            end do
        end do
    end do
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                write(12,*) v(i, j, k)!y方向速度
            end do
        end do
    end do
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                write(12,*) w(i, j, k)!z方向速度
            end do
        end do
    end do
    do k = 1, kf
        do j = 1, jf
            do i = 1, if
                write(12,*) T(i, j, k)!温度
            end do
        end do
    end do
    close(12)
end subroutine graphic