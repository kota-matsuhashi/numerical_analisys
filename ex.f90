program power_method_with_rayleigh
    implicit none
    integer, parameter :: n = 3
    real(8), dimension(n, n) :: a
    real(8), dimension(n) :: x, x_new
    real(8) :: lambda, lambda_old, tol
    integer :: max_iter, iter

    ! �p�����[�^�ݒ�
    tol = 1.0e-8
    max_iter = 1000

    ! �s��̏�����
    a = reshape((/ 2.0d0, 2.0d0, -1.0d0, &
                 0.0d0, -1.0d0, 0.0d0, &
                 0.0d0, -5.0d0, 3.0d0 /), shape(a))

    ! �����x�N�g���̐ݒ�
    x = (/ 1.0d0, 1.0d0, 1.0d0 /)
    lambda_old = 1.0d0

    do iter = 1, max_iter
        ! �s��ƃx�N�g���̐� A * x �̌v�Z�i���K���O�̃x�N�g�����g�p�j
        x_new = matmul(a, x)

        ! Rayleigh�����v�Z
        lambda = dot_product(x, x_new) / dot_product(x, x)

        ! ��������
        if (abs(lambda - lambda_old) < tol) exit
        lambda_old = lambda

        ! �ŗL�x�N�g���̍X�V�i�ő听���Ő��K���j
        x = x_new / x_new(1)  ! �ő听���Ő��K��

        ! ���ʂ̏o��
        print *, "Iteration:", iter
        print *, "Estimated eigenvalue:", lambda
        print *, "Estimated eigenvector:", x
    end do

    ! �ŏI���ʂ̏o��
    print *, "Converged eigenvalue:", lambda
    print *, "Corresponding eigenvector:", x

end program power_method_with_rayleigh
