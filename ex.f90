program power_method_with_rayleigh
    implicit none
    integer, parameter :: n = 3
    real(8), dimension(n, n) :: a
    real(8), dimension(n) :: x, x_new
    real(8) :: lambda, lambda_old, tol
    integer :: max_iter, iter

    ! パラメータ設定
    tol = 1.0e-8
    max_iter = 1000

    ! 行列の初期化
    a = reshape((/ 2.0d0, 2.0d0, -1.0d0, &
                 0.0d0, -1.0d0, 0.0d0, &
                 0.0d0, -5.0d0, 3.0d0 /), shape(a))

    ! 初期ベクトルの設定
    x = (/ 1.0d0, 1.0d0, 1.0d0 /)
    lambda_old = 1.0d0

    do iter = 1, max_iter
        ! 行列とベクトルの積 A * x の計算（正規化前のベクトルを使用）
        x_new = matmul(a, x)

        ! Rayleigh商を計算
        lambda = dot_product(x, x_new) / dot_product(x, x)

        ! 収束判定
        if (abs(lambda - lambda_old) < tol) exit
        lambda_old = lambda

        ! 固有ベクトルの更新（最大成分で正規化）
        x = x_new / x_new(1)  ! 最大成分で正規化

        ! 結果の出力
        print *, "Iteration:", iter
        print *, "Estimated eigenvalue:", lambda
        print *, "Estimated eigenvector:", x
    end do

    ! 最終結果の出力
    print *, "Converged eigenvalue:", lambda
    print *, "Corresponding eigenvector:", x

end program power_method_with_rayleigh
