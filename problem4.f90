program power_method
    implicit none
    integer, parameter :: n = 3             ! �s��̃T�C�Y
    real(8), dimension(n, n) :: a           ! ���͍s��
    real(8), dimension(n) :: x, x_new       ! �ŗL�x�N�g��
    real(8) :: lambda, lambda_old, tol      ! �ŗL�l�Ƌ��e�덷
    integer :: max_iter, iter, i
    real(8) :: norm                         ! ���[�N���b�h�m����

    ! �p�����[�^�ݒ�
    tol = 1.0e-8       ! ���e�덷�i����������j
    max_iter = 2000    ! �ő唽���񐔂𑝂₷

    ! �s��̏�����
    a = reshape((/ 2.0d0, 2.0d0, -1.0d0, & 
                 0.0d0,  -1.0d0, 0.0d0, & 
                 0.0d0, -5.0d0,  3.0d0 /), shape(a))

    ! �����x�N�g���̐ݒ�
    x = (/ 1.0d0, 1.0d0, 1.0d0 /)

    lambda = 0.0d0
    lambda_old = 1.0d0

    do iter = 1, max_iter
        ! �s��ƃx�N�g���̐� A * x �̌v�Z
        x_new = matmul(a, x)

        ! ��ꐬ���Ŋ����Đ��K���i��ꐬ�����[���łȂ����Ƃ��m�F�j
        if (x_new(1) /= 0.0d0) then
            x_new = x_new / x_new(1)
        else
            print *, "��ꐬ�����[���ł��B���K���Ɏ��s���܂����B"
            stop
        end if

        ! �ŗL�l�̐��� (Rayleigh�����g�p)
        lambda = dot_product(x_new, matmul(a, x_new)) / dot_product(x_new, x_new)

        ! ��������i���e�덷�������������ݒ�j
        if (abs(lambda - lambda_old) < tol) exit
        lambda_old = lambda

        ! �ŗL�x�N�g�� x �̍X�V
        x = x_new

        ! �e�����ł̌ŗL�l�ƌŗL�x�N�g���̏o��
        print *, "������: ", iter
        print "(A, F10.3)", "���肳�ꂽ�ŗL�l: ", lambda
        print *, "���肳�ꂽ�ŗL�x�N�g��:"
        print "(3F10.3)", x
    end do

    ! �ŏI�I�Ȍ��ʂ̏o��
    print *, "���������ŗL�l:"
    print "(F10.3)", lambda
    print *, "�Ή�����ŗL�x�N�g��:"
    print "(3F10.3)", x

end program power_method
