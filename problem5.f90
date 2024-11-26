program jacobi_method
    implicit none
    integer, parameter :: n = 4             ! 行列のサイズ
    real(8), dimension(n, n) :: a           ! 入力の対称行列
    real(8), dimension(n, n) :: eigenvectors ! 固有ベクトル行列
    real(8), dimension(n) :: eigenvalues     ! 固有値
    integer :: i, j, max_iter
    real(8) :: tol

    ! パラメータ設定
    tol = 1.0e-10      ! 許容誤差
    max_iter = 100    ! 最大反復回数

    ! 対称行列の初期化
    a = reshape([1.0d0, 1.0d0, 4.0d0, 0.0d0, &
                1.0d0, 5.0d0, 1.0d0, 0.0d0, &
                4.0d0, 1.0d0, 1.0d0, 0.0d0, &
                0.0d0, 0.0d0, 0.0d0, 2.0d0], shape(a))

    ! 固有ベクトルの初期化
    eigenvectors = 0.0d0
    do i = 1, n
        eigenvectors(i, i) = 1.0d0
    end do

    ! ヤコビ法で固有値・固有ベクトルを計算
    call jacobi(n, a, eigenvalues, eigenvectors, tol, max_iter)

    ! 最終結果出力
    print *, "最終固有値:"
    do i = 1, n
        print *, eigenvalues(i)
    end do

    print *, "最終固有ベクトル:"
    do i = 1, n
        print "(3f10.6)", eigenvectors(i, :)
    end do

contains

    ! ヤコビ法で対称行列の固有値・固有ベクトルを計算するサブルーチン
    subroutine jacobi(n, a, eigenvalues, eigenvectors, tol, max_iter)
        integer, intent(in) :: n, max_iter
        real(8), intent(in) :: tol
        real(8), dimension(n, n), intent(inout) :: a
        real(8), dimension(n), intent(out) :: eigenvalues
        real(8), dimension(n, n), intent(inout) :: eigenvectors
        integer :: iter, p, q, i
        real(8) :: theta, t, c, s, app, aqq, apq

        do iter = 1, max_iter
            ! 最大の非対角成分を探索
            call max_off_diagonal(n, a, p, q)
            if (abs(a(p, q)) < tol) exit

            ! 回転角度の計算
            theta = 0.5d0 * atan2(2.0d0 * a(p, q), a(q, q) - a(p, p))
            c = cos(theta)
            s = sin(theta)

            ! 行列Aの更新
            app = a(p, p)
            aqq = a(q, q)
            apq = a(p, q)
            a(p, p) = c**2 * app - 2.0d0 * s * c * apq + s**2 * aqq
            a(q, q) = s**2 * app + 2.0d0 * s * c * apq + c**2 * aqq
            a(p, q) = 0.0d0
            a(q, p) = 0.0d0

            do i = 1, n
                if (i /= p .and. i /= q) then
                    apq = a(i, p)
                    a(i, p) = c * apq - s * a(i, q)
                    a(p, i) = a(i, p)
                    a(i, q) = s * apq + c * a(i, q)
                    a(q, i) = a(i, q)
                end if
            end do

            ! 固有ベクトルの更新
            do i = 1, n
                apq = eigenvectors(i, p)
                eigenvectors(i, p) = c * apq - s * eigenvectors(i, q)
                eigenvectors(i, q) = s * apq + c * eigenvectors(i, q)
            end do

            ! 固有値の更新（対角成分）
            do i = 1, n
                eigenvalues(i) = a(i, i)
            end do

            ! 各反復ステップで固有値を出力
            print *, "反復回数:", iter
            print *, "固有値:"
            do i = 1, n
                print *, eigenvalues(i)
            end do
        end do
    end subroutine jacobi

    ! 非対角成分の最大値を探すサブルーチン
    subroutine max_off_diagonal(n, a, p, q)
        integer, intent(in) :: n
        real(8), dimension(n, n), intent(in) :: a
        integer, intent(out) :: p, q
        integer :: i, j
        real(8) :: max_val

        max_val = 0.0d0
        p = 1
        q = 2
        do i = 1, n-1
            do j = i+1, n
                if (abs(a(i, j)) > max_val) then
                    max_val = abs(a(i, j))
                    p = i
                    q = j
                end if
            end do
        end do
    end subroutine max_off_diagonal

end program jacobi_method
