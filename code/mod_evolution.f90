module evolution

  use pf_setting
  implicit none

contains

  SUBROUTINE init_g_T0()
    implicit none
    ! calculate without B_string tech from time = 1 to ntime, yielding green function at time = ntime+1
    integer :: i, p, block_num, i_pf , i_pla
    complex(8), allocatable :: lmat(:, :), rmat(:, :), PBP(:, :)
    complex(8),allocatable :: lmat_T(:,:), rmat_T(:, :), PBP_T(:, :)
    type(pf_data_type),pointer :: p_data
    IF (.not. allocated(g)) allocate (g(Ns, Ns))
    IF (.not. allocated(g_h)) allocate (g_h(Ns, Ns))
    IF (.not. allocated(Q_string)) allocate (Q_string(Ns, nelec, nblock),&
    & D_string(nelec, nblock))

    !print*,'slater_D:',slater_D
    ln_cw = 0d0 ! initialize configuration weight
    allocate (rmat(Ns, nelec), lmat(nelec, Ns), PBP(nelec, nelec))

         !! initialize QDR_string(block = 1 and nblock)

    Q_string(:, 1:nelec, 1) = slater_Q(:, 1:nelec)
    Q_string(:, 1:nelec, nblock) = slater_Q(:, 1:nelec)
    if(second_trotter) then
      !Q_string(:, 1:nelec, 1) =  matmul(expK_half , Q_string(:, 1:nelec, 1))
      !Q_string(:, 1:nelec, nblock) =  matmul(expK_inv_half,Q_string(:, 1:nelec, nblock))
    end if
    call qdr(ns, nelec, Q_string(:, 1:nelec, 1),  &
    & R_string(1:nelec, 1:nelec),D_string(1:nelec, 1))
    call qdr(ns, nelec, Q_string(:, 1:nelec, nblock),  &
   & R_string(1:nelec, 1:nelec),D_string(1:nelec, nblock))
    rmat = Q_string(:, 1:nelec, 1)

    DO p = 1, ntime
      do i_pf = 1, n_phonon_field
        ! left evolve the rmat matrix : expK * expV * Q
        call left_evolve(p, rmat, nelec,i_pf,two_way=.false.)

      end do

      if (mod(p, ngroup) /= 0) cycle; 
      block_num = (p - 1)/ngroup + 2
      do i = 1, nelec

        !rmat(:,i) = rmat(:,i) * D_string(i,block_num-1)
      end do
      call qdr(ns, nelec, rmat, R_string(1:nelec, 1:nelec), D_string(1:nelec, block_num))

      Q_string(:, 1:nelec, block_num) = rmat
      !print*,'D_string:',block_num,':',D_string(1:nelec,flv,block_num)
      !R_string(1:nelec,1:nelec,flv,block_num) = matmul(R_string(1:nelec,1:nelec,flv,block_num),&
      !& R_string(1:nelec,1:nelec,flv,block_num-1))

    end do

    !ln_cw = ln_cw + sum(log(D_string(1:nelec,flv,nblock-1)))
         !! abs only for sign problem free to avoid minus value of D_string

    lmat = conjg(transpose(Q_string(:, 1:nelec, nblock)))
    call cal_g_and_det(lmat,rmat)

    

    deallocate (rmat, lmat, PBP)

    !print*,'max g:',maxval(abs(g_h))

  end SUBROUTINE

  SUBROUTINE get_g_scratch_T0(time, forward, loop)
    implicit none
    ! come here when mod(time_count,ngroup)) == 0 including time = ntime
    ! when forward, should give g(time+1)
    ! when backward, should give g(time)
    integer, intent(in) :: time, loop
    logical, intent(in) :: forward
    integer :: i, p, block_num,i_pf
    complex(8), allocatable :: lmat(:, :), rmat(:, :), PBP(:, :)
    complex(8) :: gfast(ns, ns)
    complex(8) :: ln_cw_fast
    ln_cw_fast = ln_cw
    ln_cw = 0d0
    !print*,'qmat:',qmat

    allocate (rmat(Ns, nelec), lmat(nelec, Ns), PBP(nelec, nelec))

    gfast = g(:, :)

    if (forward) then

      block_num = time/ngroup + 1
      rmat = Q_string(:, 1:nelec, block_num - 1)

      DO p = time - ngroup, time - 1
        do i_pf = 1 , n_phonon_field
          call left_evolve(p, rmat, nelec,i_pf,two_way=.false.) !expK * expV * Q
        end do
      end do
      !forall(i = 1:nelec) rmat(:,i) = rmat(:,i) * D_string(i,flv,block_num-1) !Q' = Q*D
      CALL qdr(ns, nelec, rmat, R_string(1:nelec, 1:nelec), D_string(1:nelec, block_num))
      Q_string(:, 1:nelec, block_num) = rmat
      !R_string(1:nelec,1:nelec,flv,block_num) = matmul(R_string(1:nelec,1:nelec,flv,block_num),&
      !& R_string(1:nelec,1:nelec,flv,block_num-1)) !
      lmat = conjg(transpose(Q_string(:, 1:nelec, block_num + 1)))
      !update ln_cw
      
      !ln_cw = ln_cw + sum(log(D_string(1:nelec,flv,block_num+1)))
    else !evolve backwards

      block_num = time/ngroup + 2  ! depends on size of ngroup
      !print*,"time,block_num:",time,block_num
      if (ngroup == 1) block_num = block_num - 1
      lmat = conjg(transpose(Q_string(:, 1:nelec, block_num + 1)))

      DO p = time + ngroup - 1, time, -1
        do i_pf = n_phonon_field,1,-1
          call right_evolve(p, lmat, nelec,i_pf,two_way=.false.)!Q*expK * expV
        end do
      end do

      ! forall(i = 1:nelec) lmat(i,:) = lmat(i,:) * D_string(i,flv,block_num+1) !Q' = Q*D

      CALL ldq(nelec, Ns, lmat, R_string(1:nelec, 1:nelec), D_string(1:nelec, block_num))

      Q_string(:, 1:nelec, block_num) = conjg(transpose(lmat))
      D_string(1:nelec, block_num) = D_string(1:nelec, block_num)
      !R_string(1:nelec,1:nelec,flv,block_num) = matmul(R_string(1:nelec,1:nelec,flv,block_num+1),&
      !& R_string(1:nelec,1:nelec,flv,block_num)) ! attention that R_string was not hermitian_conjugated
      rmat = Q_string(:, 1:nelec, block_num - 1)
      !ln_cw = ln_cw + sum(log(D_string(1:nelec,flv,block_num-1)))
    end if
    ! calculate g =  1 - rmat * (lmat*rmat)^(-1) * lmat
    call cal_g_and_det(lmat,rmat)
    ! calculate difference between fast updating g and scratch
    ! difference
    err_fast = maxval(abs(gfast - g(:, :)))
    if (err_fast > err_fast_max) then
      print *, 'time:', time, 'error_fast: ', err_fast, 'loop:', loop, 'myid:', myid
      err_fast_max = err_fast
      print *, 'gfast:', maxval(abs(gfast))
      print *, 'g:', maxval(abs(g(:, :)))
      !stop 'Error in update: difference too large'
    end if
    err_fast = abs((real(ln_cw - ln_cw_fast)))
    if (err_fast > err_fast_max) then
      err_fast_max = err_fast
      print *, 'time:', time, 'error_fast_det: ', err_fast, 'loop:', loop, 'myid:', myid
      print *, 'lnw_fast,lnw_from_scratch:', ln_cw_fast, ln_cw
      !stop 'Error in update lnw: difference too large'
    end if
    deallocate (rmat, lmat, PBP)
  end SUBROUTINE

  subroutine get_g_scratch_tau(time,rmat,mat)
    implicit none
    !> calculate gf from scratch especially for tau and update the rmat matrix
    ! like get_g_scratch_T0, but only for tau
    ! input time+1 to give g_h(time+1)
    ! e.g. time = 46
    integer, intent(in) :: time
    complex(8), intent(inout) :: mat(Ns, Ns)
    complex(8) :: gfast(Ns, Ns)
    complex(8) :: PBP(nelec, nelec)
    complex(8), intent(inout) :: rmat(Ns, nelec)
    complex(8) :: lmat(nelec, Ns), D_string_temp(nelec)
    integer :: i, p, block_num,i_pf
    block_num = time/ngroup + 1 ! e.g. block_num = 46/5 + 1 = 10

    gfast = mat
    DO p = time - ngroup, time - 1
      do i_pf =  n_phonon_field,1,-1
        call left_evolve(p, rmat, nelec,i_pf,two_way=.false.) !expK *  expV * Q
      end do
    end do
    CALL qdr(ns, nelec, rmat, R_string(1:nelec, 1:nelec), D_string_temp)
    lmat = conjg(transpose(Q_string(:, 1:nelec, block_num + 1))) ! e.g. block_num + 1 = 11
    PBP = matmul(lmat, rmat) !lmat*rmat->PBP
    call inverse(nelec, PBP)
    mat = MATMUL(rmat, matmul(PBP, lmat))

    err_fast = maxval(abs(gfast - mat))

    if (err_fast > err_fast_max) then

      print*, 'fast_update_error in get_g_scratch_tau:'
      print *, 'time:', time, 'error_fast: ', err_fast,  'myid:', myid
      err_fast_max = err_fast
      !stop 'Error in tau-update: difference too large'
    end if

  end subroutine get_g_scratch_tau

  SUBROUTINE left_evolve(time, mat, n,i_pf,two_way)
    !mat -> expK1*expK2...expKn*mat*exp-Kn...*exp-K2*exp-K1
    IMPLICIT none
    integer, intent(in) :: time, n,i_pf
    logical, intent(in):: two_way
    complex(8) :: mat(ns, n)
    integer :: i_pla
    type(pf_type),pointer :: ppf
    type(pf_data_type), pointer :: p_data
    ppf => pf_list(i_pf)
    p_data => ppf%p_data
    do i_pla = 1, ppf%n_plaquette
      mat(p_data%pla_site_list(:,i_pla), :) = matmul(p_data%expKV(:,:,i_pla,time) , mat(p_data%pla_site_list(:,i_pla), :))
      if (two_way) then
        if (n /= ns) stop 'error in left_evolve, dimension not for evolve two-way'
        mat(:, p_data%pla_site_list(:,i_pla)) = matmul(mat(:, p_data%pla_site_list(:,i_pla)), p_data%expKV_inv(:,:,i_pla,time))
      end if
    end do

  end SUBROUTINE left_evolve

  SUBROUTINE right_evolve(time, mat, n, i_pf,two_way)
    IMPLICIT none
    integer, intent(in) :: time, n,i_pf
    logical, intent(in) :: two_way
    complex(8) :: mat(n,ns)
    integer :: i_pla
    type(pf_type),pointer :: ppf
    type(pf_data_type), pointer :: p_data
    ppf => pf_list(i_pf)
    p_data => ppf%p_data
    do i_pla = 1, ppf%n_plaquette
      mat(:,p_data%pla_site_list(:,i_pla)) = matmul(mat(:,p_data%pla_site_list(:,i_pla)),p_data%expKV(:,:,i_pla,time))
      if (two_way) then
        if (n /= ns) stop 'error in right_evolve, dimension not for evolve two-way'
        mat(p_data%pla_site_list(:,i_pla), :) = matmul(p_data%expKV_inv(:,:,i_pla,time),mat(p_data%pla_site_list(:,i_pla), :))
      end if
    end do

  end SUBROUTINE right_evolve

  subroutine cal_g_and_det(lmat,rmat)
    implicit none
    integer :: i
    complex(8), intent(in) :: lmat(nelec,Ns), rmat(Ns,nelec)
    complex(8),allocatable :: PBP(:,:)
    complex(8) :: lmat_T(2*nelec,Ns), rmat_T(Ns,2*nelec)
    complex(8) :: inner_prod(nelec,nelec)
    ln_cw = 0d0
    if(.true.) then
      allocate(PBP(nelec,nelec))
      PBP = matmul(lmat,rmat)
      ln_cw = ln_cw + log((det(nelec, PBP))) !> abs for real only
      ln_cw = ln_cw + sum(log((D_string)))
      call inverse(nelec, PBP) ! PBP^(-1) -> PBP
      g_h(:, :) = MATMUL(rmat, matmul(PBP, lmat))
    else
      allocate(PBP(2*nelec,2*nelec))
      rmat_T(:,1:nelec) = rmat(:,1:nelec)
      rmat_T(:,nelec+1:2*nelec) = conjg(MATMUL(TR_mat,rmat(:,1:nelec)))
      lmat_T(1:nelec,:) = lmat(1:nelec,:)
      lmat_T(nelec+1:2*nelec,:) = conjg(MATMUL(lmat(1:nelec,:),transpose(TR_mat)))
      PBP = matmul(lmat_T,rmat_T)
      do i = 1, 2*nelec
        !print*,'rmat:',i,rmat_T(:,i)
      end do
      inner_prod = matmul(conjg(transpose(Rmat_T(:,1: nelec))),lmat_T(:,1: nelec))
      do i = 1 , nelec
         print*,'max inner_prod:', abs(inner_prod(i,:))
      end do
      ! print*,'inner_prod:',inner_prod
      print*,'det inner_prod:',det(nelec,inner_prod)
      print*,'det PBP:',det(2*nelec,PBP)
      
      ln_cw = ln_cw + log(abs(det(2*nelec, PBP))) !> abs for real only      
      ln_cw = ln_cw + sum(log(abs(D_string))) * 2
      call inverse(2*nelec, PBP) ! PBP^(-1) -> PBP
      g_h(:, :) = MATMUL(rmat_T, matmul(PBP, lmat_T))
      stop
    end if
    g(:, :) = -g_h(:, :)
    forall (i=1:Ns) g(i, i) = g(i, i) + 1d0
    ln_cw = ln_cw*ncopy
    ln_cw = complex(real(ln_cw),mod(aimag(ln_cw),2*pi)) ! make sure weight is in the principal branch
  end subroutine cal_g_and_det
end module evolution