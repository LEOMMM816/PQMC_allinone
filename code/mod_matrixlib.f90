MODULE matrixlib
  
  INTERFACE trace
    MODULE PROCEDURE dtrace, ztrace
  END INTERFACE

  INTERFACE qdr ! QDR' decomposition from QR decomposition by DR' = R, R'(i,i) = 1
    MODULE PROCEDURE dqdr, zqdr
  END INTERFACE

  INTERFACE ldq
    MODULE PROCEDURE dldq, zldq
  END INTERFACE

  INTERFACE eigen
    MODULE PROCEDURE deigen, zeigen
  END INTERFACE

  INTERFACE inverse
    MODULE PROCEDURE dinverse, zinverse
  END INTERFACE

  INTERFACE qr !QR decomposition
    MODULE PROCEDURE dqr, zqr
  END INTERFACE

  INTERFACE lq !LQ decomposition
    MODULE PROCEDURE dlq, zlq
  END INTERFACE

  INTERFACE det
    MODULE PROCEDURE ddet, zdet
  END INTERFACE

  INTERFACE expm
    MODULE PROCEDURE dexpm, zexpm
  END INTERFACE

  INTERFACE sort
    MODULE PROCEDURE rsort, dsort
  END INTERFACE

  INTERFACE swap
    MODULE PROCEDURE rswap, dswap, cswap
  END INTERFACE
CONTAINS
!--------------------------------------------------
!   general functins
!--------------------------------------------------

  integer function modi(n, p)
    integer, intent(in)::n, p

    modi = mod(n, p)
    if (n < 0) then
      modi = p + modi
      return

    else if (modi == 0) then
      modi = p
    end if
  end function

  real(8) function rlog(x) result(re)
    real(8) :: x

    re = log(abs(x))

  end function

  function idx(n) result(re)
    integer, intent(in) :: n
    integer :: i
    real(8) :: re(n, n)
    re = 0
    forall (i=1:n) re(i, i) = 1d0
    return
  end function

  function idx_c(n) result(re)
    integer, intent(in) :: n
    integer :: i
    complex(8) :: re(n, n)
    re = 0
    forall (i=1:n) re(i, i) = (1d0, 0d0)
    return
  end function

  function sigma_x(n) result(re)
    integer, intent(in) :: n

    real(8) :: re(n, n)
    re = 0
    if (n /= 2) stop 'n/=2 in sigma_x'
    re(1, 2) = 1d0
    re(2, 1) = 1d0
    return
  end function

  function sigma_y(n) result(re) ! actually i*sigma_y
    integer, intent(in) :: n

    complex(8) :: re(n, n)
    re = 0
    if (n /= 2) stop 'n/=2 in sigma_x'
    re(1, 2) = 1
    re(2, 1) = -1
    return
  end function

  function sigma_z(n) result(re)
    integer, intent(in) :: n

    real(8) :: re(n, n)
    re = 0
    if (n /= 2) stop 'n/=2 in sigma_x'
    re(1, 1) = 1d0
    re(2, 2) = -1d0
    return
  end function
  subroutine rsort(array_, n, ascending)
    implicit none
    integer, intent(in) :: n
    real(8) :: array_(n)
    logical :: ascending
    integer :: i, ind

    if (ascending) then
      do i = 1, n
        ind = minloc(array_(i:n), dim=1) + i - 1
        if (array_(ind) /= array_(i)) then
          call swap(array_(ind), array_(i))
        end if
      end do
    else
      do i = 1, n
        ind = maxloc(array_(i:n), dim=1) + i - 1
        if (array_(ind) /= array_(i)) then
          call swap(array_(ind), array_(i))
        end if
      end do
    end if
  end subroutine

  subroutine dsort(array_, n, ascending)
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: array_(n)
    logical :: ascending
    integer :: i, ind

    if (ascending) then
      do i = 1, n
        ind = minloc(array_(i:n), dim=1) + i - 1
        if (array_(ind) /= array_(i)) then
          call swap(array_(ind), array_(i))
        end if
      end do
    else
      do i = 1, n
        ind = maxloc(array_(i:n), dim=1) + i - 1
        if (array_(ind) /= array_(i)) then
          call swap(array_(ind), array_(i))
        end if
      end do
    end if
  end subroutine
  subroutine rswap(a, b)
    implicit none
    real(8) :: a, b, temp
    temp = a
    a = b
    b = temp
  end subroutine

  subroutine dswap(a, b)
    implicit none
    integer :: a, b, temp
    temp = a
    a = b
    b = temp
  end subroutine

  subroutine cswap(a, b)
    implicit none
    complex(8) :: a, b, temp
    temp = a
    a = b
    b = temp
  end subroutine

  subroutine dexpm(mat, n)
    implicit none
    integer, intent(in) :: n
    logical :: Hermitian
    real(8) :: v(n)
    real(8) :: mat(n, n), exp_mat(n, n)
    complex(8) :: c_mat(n, n)
    integer :: i, j
    if(n==0) stop 'n=0 in dexpm'
    if(n==1) then
      mat(1,1) = exp(mat(1,1))
      return
    end if
    c_mat = transpose(mat)
    if(maxval(abs(mat - c_mat)) < 1d-7) then
      Hermitian = .true.
    else if( maxval(abs(mat + c_mat)) < 1d-7) then
      Hermitian = .false.
    else
      stop 'matrix is neither symmetric nor anti-symmetric in dexpm'
    end if
    if(Hermitian) then
      call eigen(n, mat, v)
      do j = 1, n
        do i = 1, n
          exp_mat(i, j) = sum(mat(i, :)*exp(v(:))*mat(j, :))
        end do
      end do
    else
      ! anti-Hermitian matrix
      c_mat = complex(0d0,1d0) * mat
      call eigen(n, c_mat, v)
      do j = 1, n
        do i = 1, n
          exp_mat(i, j) = real(sum(c_mat(i, :)*exp(complex(0d0,-1d0) * v(:))*conjg(c_mat(j, :))))
        end do
      end do

    end if
    mat = exp_mat
    return
  end subroutine
  subroutine zexpm(mat, n)
    implicit none
    integer, intent(in) :: n
    logical :: Hermitian
    real(8) :: v(n)
    complex(8) :: mat(n, n), exp_mat(n, n),c_mat(n,n)
    integer :: i, j
    if(n==0) stop 'n=0 in dexpm'
    if(n==1) then
      mat(1,1) = exp(mat(1,1))
      return
    end if
    c_mat = transpose(conjg(mat))
    if(maxval(abs(mat - c_mat)) < 1d-7) then
      Hermitian = .true.
    else if(maxval(abs(mat + c_mat)) < 1d-7) then
      Hermitian = .false.
    else
      stop 'matrix is neither Hermitian nor anti-Hermitian in zexpm'
    end if
    if(Hermitian) then
      call eigen(n, mat, v)
      do j = 1, n
        do i = 1, n
          exp_mat(i, j) = sum(mat(i, :)*exp(v(:))*conjg(mat(j, :)))
        end do
      end do
    else
      ! anti-Hermitian matrix
      c_mat = complex(0d0,1d0) * mat
      call eigen(n, c_mat, v)
      do j = 1, n
        do i = 1, n
          exp_mat(i, j) = sum(mat(i, :)*exp(complex(0d0,-1d0) * v(:))*conjg(mat(j, :)))
        end do
      end do
    end if
    mat = exp_mat
  end subroutine
  FUNCTION dtrace(n, a)
    IMPLICIT NONE
    INTEGER n, i
    REAL(8) a(n, n), dtrace
    dtrace = 0d0
    DO i = 1, n
      dtrace = dtrace + a(i, i)
    END DO
  END FUNCTION

  FUNCTION ztrace(n, a)
    IMPLICIT NONE
    INTEGER n, i
    COMPLEX(8) a(n, n), ztrace
    ztrace = (0d0, 0d0)
    DO i = 1, n
      ztrace = ztrace + a(i, i)
    END DO
  END FUNCTION

  SUBROUTINE dqdr(m, n, a, r, d)
    IMPLICIT NONE
    INTEGER m, n, i
    REAL(8) a(m, n), r(n, n), d(n)

    CALL dqr(m, n, a, r)
    DO i = 1, n
      d(i) = r(i, i)
      r(i, i:n) = r(i, i:n)/d(i)
    END DO

  END SUBROUTINE

  SUBROUTINE zqdr(m, n, a, r, d)
    IMPLICIT NONE
    INTEGER m, n, i
    COMPLEX(8) a(m, n), r(n, n), d(n)

    CALL zqr(m, n, a, r)
    DO i = 1, n
      d(i) = r(i, i)
      r(i, i:n) = r(i, i:n)/d(i)
    END DO

  END SUBROUTINE

  SUBROUTINE dldq(m, n, a, l, d)
    IMPLICIT NONE
    INTEGER m, n, i
    REAL(8) a(m, n), l(m, m), d(m)

    CALL dlq(m, n, a, l)
    DO i = 1, m
      d(i) = l(i, i)
      l(i:m, i) = l(i:m, i)/d(i)
    END DO

  END SUBROUTINE

  SUBROUTINE zldq(m, n, a, l, d)
    IMPLICIT NONE
    INTEGER m, n, i
    COMPLEX(8) a(m, n), l(m, m), d(m)

    CALL zlq(m, n, a, l)
    DO i = 1, m
      d(i) = l(i, i)
      l(i:m, i) = l(i:m, i)/d(i)
    END DO

  END SUBROUTINE

  SUBROUTINE deigen(n, a, v)
    IMPLICIT NONE
    INTEGER n, info
    REAL(8) a(n, n), v(n), work(64*n)
    IF (n <= 0) STOP 'N<=0 is invalid, inside DEIGEN'
    IF (n == 1) THEN
      v(1) = a(1, 1); a(1, 1) = 1d0; RETURN
    END IF
    CALL dsyev('V', 'U', n, a, n, v, work, 64*n, info)
    IF (info /= 0) STOP 'ERROR @ DSYEV, inside DEIGEN'
  END SUBROUTINE

  SUBROUTINE zeigen(n, a, v)
    IMPLICIT NONE
    INTEGER n, info
    REAL(8) rwork(3*n - 2), v(n)
    COMPLEX(8) work(64*n), a(n, n)
    IF (n <= 0) STOP 'N<=0 is invalid, inside DEIGEN'
    IF (n == 1) THEN
      v(1) = real(a(1, 1)); a(1, 1) = 1d0; RETURN
    END IF
    CALL zheev('v', 'u', n, a, n, v, work, 64*n, rwork, info)
    IF (info /= 0) STOP 'ERROR @ ZHEEV, inside ZEIGEN'
  END SUBROUTINE

  SUBROUTINE dinverse(n, a)
    IMPLICIT NONE
    INTEGER n, info, ipiv(n)
    REAL(8) a(n, n), work(64*n), b(2, 2)
    IF (n <= 0) STOP 'N<=0 is invalid, inside DINVERSE'
    IF (n == 1) THEN
      a(1, 1) = 1d0/a(1, 1); RETURN
    END IF
    !IF(n==2)THEN
    !   b(1,1) = a(2,2)
    !   b(2,2) = a(1,1)
    !   b(1,2) = -a(2,1)
    !   b(2,1) = -a(1,2)
    !   b = b / det(n,a)
    !   a = b
    !   return
    !END IF
    CALL dgetrf(n, n, a, n, ipiv, info)
    IF (info /= 0) STOP 'ERROR @ DGETRF, inside DINVERSE'
    CALL dgetri(n, a, n, ipiv, work, 64*n, info)
    IF (info /= 0) STOP 'ERROR @ DGETRI, inside DINVERSE'
  END SUBROUTINE

  SUBROUTINE zinverse(n, a)
    IMPLICIT NONE
    INTEGER n, info, ipiv(n)
    COMPLEX(8) a(n, n), work(64*n), b(2, 2)
    IF (n <= 0) STOP 'N<=0 is invalid, inside ZINVERSE'
    IF (n == 1) THEN
      a(1, 1) = 1d0/a(1, 1); RETURN
    END IF
    !IF(n==2)THEN
    !   b(1,1) = a(2,2)
    !   b(2,2) = a(1,1)
    !   b(1,2) = -a(2,1)
    !   b(2,1) = -a(1,2)
    !   b = b / det(n,a)
    !   a = b
    !   return
    !END IF
    CALL zgetrf(n, n, a, n, ipiv, info)
    IF (info /= 0) STOP 'ERROR @ ZGETRF, inside ZINVERSE'
    CALL zgetri(n, a, n, ipiv, work, 64*n, info)
    IF (info /= 0) STOP 'ERROR @ ZGETRI, inside ZINVERSE'
  END SUBROUTINE

  SUBROUTINE dqr(m, n, a, r)
    IMPLICIT NONE
    INTEGER m, n, info, j, lwork
    REAL(8) a(m, n), r(n, n), tau(n), work(n*64)

    IF (m < n) STOP 'm<n is invalid, inside DQR'
    IF (n <= 0) STOP 'n<=0 is invalid, inside DQR'
    lwork = n*64
    CALL dgeqrf(m, n, a, m, tau, work, lwork, info)
    IF (info /= 0) STOP 'ERROR @ ZGEQRF, inside DQR'
    r = 0d0
    DO j = 1, n
      r(1:j, j) = a(1:j, j)
    END DO
    CALL dorgqr(m, n, n, a, m, tau, work, lwork, info)
    IF (info /= 0) STOP 'ERROR @ ZUNGQR, inside DQR'

  END SUBROUTINE

  SUBROUTINE zqr(m, n, a, r)
    IMPLICIT NONE
    INTEGER m, n, info, j, lwork
    COMPLEX(8) a(m, n), r(n, n), tau(n), work(n*64)
    IF (m < n) STOP 'm<n is invalid, inside ZQR'
    IF (n <= 0) STOP 'n<=0 is invalid, inside ZQR'
    lwork = n*64
    CALL zgeqrf(m, n, a, m, tau, work, lwork, info)
    IF (info /= 0) STOP 'ERROR @ ZGEQRF, inside ZQR'
    r = 0d0
    DO j = 1, n
      r(1:j, j) = a(1:j, j)
    END DO
    CALL zungqr(m, n, n, a, m, tau, work, lwork, info)
    IF (info /= 0) STOP 'ERROR @ ZUNGQR, inside ZQR'
  END SUBROUTINE

  SUBROUTINE dlq(m, n, a, l)
    IMPLICIT NONE
    INTEGER m, n, info, j, lwork
    REAL(8) a(m, n), l(m, m), tau(m), work(m*64)

    IF (m > n) STOP 'm>n is invalid, inside DLQ'
    IF (m <= 0) STOP 'n<=0 is invalid, inside DLQ'
    lwork = m*64
    CALL dgelqf(m, n, a, m, tau, work, lwork, info)
    IF (info /= 0) STOP 'ERROR @ ZGELQF, inside DLQ'
    l = 0d0
    DO j = 1, m
      l(j:m, j) = a(j:m, j)
    END DO
    CALL dorglq(m, n, m, a, m, tau, work, lwork, info)
    IF (info /= 0) STOP 'ERROR @ ZUNGLQ, inside DLQ'

  END SUBROUTINE

  SUBROUTINE zlq(m, n, a, l)
    IMPLICIT NONE
    INTEGER m, n, info, j, lwork
    COMPLEX(8) a(m, n), l(m, m), tau(m), work(m*64)
    IF (m > n) STOP 'm>n is invalid, inside ZLQ'
    IF (m <= 0) STOP 'm<=0 is invalid, inside ZLQ'
    lwork = m*64
    CALL zgelqf(m, n, a, m, tau, work, lwork, info)
    IF (info /= 0) STOP 'ERROR @ ZGELQF, inside ZLQ'
    l = 0d0
    DO j = 1, m
      l(j:m, j) = a(j:m, j)
    END DO
    CALL zunglq(m, n, m, a, m, tau, work, lwork, info)
    IF (info /= 0) STOP 'ERROR @ ZUNGLQ, inside ZLQ'
  END SUBROUTINE

!========================================================================================
! IMPORTANT: there may be unstability problems in these subroutines solving determinant.
! In fact, calculating determinant is not recommendated in real calculations.
!========================================================================================

  FUNCTION ddet(n, a)
    IMPLICIT NONE
    INTEGER n, i, info, ipvt(n)
    REAL(8) a(n, n), b(n, n), ddet
    IF (n <= 0) STOP 'N<0 is invalid, inside DDET'
    IF (n == 1) THEN
      ddet = a(1, 1); RETURN
    END IF
    !IF(n==2)THEN
    !   ddet=a(1,1)*a(2,2)-a(1,2)*a(2,1);RETURN
    !END IF
    b = a
    CALL dgetrf(n, n, b, n, ipvt, info)
    IF (info /= 0) THEN
      ddet = 0d0; RETURN
    END IF
    info = 1
    ddet = 0d0
    DO i = 1, n
      IF (ipvt(i) /= i) info = -info
      IF (b(i, i) < 0d0) THEN
        info = -info
        b(i, i) = -b(i, i)
      END IF
      ddet = ddet + log(b(i, i))
    END DO
    ddet = exp(ddet)
    IF (info < 0) ddet = -ddet
  END FUNCTION

  FUNCTION zdet(n, a)
    use, intrinsic :: ieee_arithmetic
    IMPLICIT NONE
    INTEGER i, n, info, ipvt(n)
    COMPLEX(8) a(n, n), b(n, n), zdet
    IF (n <= 0) STOP 'N<0 is invalid, inside ZDET'
    IF (n == 1) THEN
      zdet = a(1, 1); RETURN
    END IF
    !IF(n==2)THEN
    !   zdet=a(1,1)*a(2,2)-a(1,2)*a(2,1);RETURN
    !END IF
    b = a
    CALL zgetrf(n, n, b, n, ipvt, info)
    IF (info /= 0) THEN
      zdet = 0d0; RETURN
    END IF
    info = 1
    zdet = 1d0
    !zdet=0d0
    DO i = 1, n
      IF (ipvt(i) /= i) info = -info
      
      zdet = zdet*b(i, i)
      !zdet=zdet+log(b(i,i))
    END DO
    !zdet=exp(zdet)
    IF (info < 0) zdet = -zdet
  END FUNCTION

END MODULE

