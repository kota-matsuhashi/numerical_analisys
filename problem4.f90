program power_method
    implicit none
    integer, parameter :: n = 3             ! 行列のサイズ
    real(8), dimension(n, n) :: a           ! 入力行列
    real(8), dimension(n) :: x, x_new       ! 固有ベクトル
    real(8) :: lambda, lambda_old, tol      ! 固有値と許容誤差
    integer :: max_iter, iter, i
    real(8) :: norm                         ! ユークリッドノルム

    ! パラメータ設定
    tol = 1.0e-8       ! 許容誤差（厳しくする）
    max_iter = 2000    ! 最大反復回数を増やす

    ! 行列の初期化
    a = reshape((/ 2.0d0, 2.0d0, -1.0d0, & 
                 0.0d0,  -1.0d0, 0.0d0, & 
                 0.0d0, -5.0d0,  3.0d0 /), shape(a))

    ! 初期ベクトルの設定
    x = (/ 1.0d0, 1.0d0, 1.0d0 /)

    lambda = 0.0d0
    lambda_old = 1.0d0

    do iter = 1, max_iter
        ! 行列とベクトルの積 A * x の計算
        x_new = matmul(a, x)

        ! 第一成分で割って正規化（第一成分がゼロでないことを確認）
        if (x_new(1) /= 0.0d0) then
            x_new = x_new / x_new(1)
        else
            print *, "第一成分がゼロです。正規化に失敗しました。"
            stop
        end if

        ! 固有値の推定 (Rayleigh商を使用)
        lambda = dot_product(x_new, matmul(a, x_new)) / dot_product(x_new, x_new)

        ! 収束判定（許容誤差を少し厳しく設定）
        if (abs(lambda - lambda_old) < tol) exit
        lambda_old = lambda

        ! 固有ベクトル x の更新
        x = x_new

        ! 各反復での固有値と固有ベクトルの出力
        print *, "反復回数: ", iter
        print "(A, F10.3)", "推定された固有値: ", lambda
        print *, "推定された固有ベクトル:"
        print "(3F10.3)", x
    end do

    ! 最終的な結果の出力
    print *, "収束した固有値:"
    print "(F10.3)", lambda
    print *, "対応する固有ベクトル:"
    print "(3F10.3)", x

end program power_method
