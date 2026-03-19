module evolution

  use pf_setting
  implicit none

  ABSTRACT INTERFACE
    SUBROUTINE init_g_f(start_time)
      integer, intent(in) :: start_time
    END SUBROUTINE init_g_f
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE get_g_scratch_f(time, forward, loop)
      integer, intent(in) :: time,loop
      logical, intent(in) :: forward
    END SUBROUTINE get_g_scratch_f
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE cal_new_weight_f()
    END SUBROUTINE cal_new_weight_f
  END INTERFACE
  ! 2. 声明过程指针

  PROCEDURE(init_g_f), POINTER :: init_g => NULL()
  procedure(get_g_scratch_f), POINTER :: get_g_scratch => NULL()
  procedure(cal_new_weight_f), POINTER :: cal_new_weight => NULL()
contains

  subroutine set_temperature()
    implicit none
    ! set PQMC or FTQMC
    if(proj) then
      ! PQMC
      init_g => init_g_T0
      get_g_scratch => get_g_scratch_T0
      cal_new_weight => cal_new_weight_T0
    else
      ! FTQMC
      init_g => init_g_FT
      get_g_scratch => get_g_scratch_FT
      cal_new_weight => cal_new_weight_FT
    end if

  end subroutine
  SUBROUTINE init_g_T0(start_time)

    implicit none
    integer, intent(in) :: start_time
    integer :: start_block, start_time_by_block
    ! recalculate B_string(s) from time = start_time to ntime, yielding green function at time = ntime+1
    ! start_time ranges from 1 to ntime, e.g. start_time = 1 means recalulate all B_strings
    ! start_time = ntime means no evolution, just calculate g_h at ntime+1
    ! start_time and start_block are where the evolution starts.
    integer ::  p, block_num, i_pf
    complex(dp), allocatable :: lmat(:, :), rmat(:, :), PBP(:, :)
    IF (.not. allocated(g)) allocate (g(Ns, Ns))
    IF (.not. allocated(g_h)) allocate (g_h(Ns, Ns))
    IF (.not. allocated(Q_string)) allocate (Q_string(Ns, nelec, nblock),D_string(nelec, nblock),R_string(nelec, nelec, nblock))

    ln_cw = 0d0 ! initialize configuration weight
    allocate (rmat(Ns, nelec), lmat(nelec, Ns), PBP(nelec, nelec))

         !! initialize QDR_string(block = 1 and nblock)

    Q_string(:, 1:nelec, 1) = slater_Q(:, 1:nelec)
    Q_string(:, 1:nelec, nblock) = slater_Q(:, 1:nelec)
    D_string(:, 1) = slater_D(1:nelec)
    D_string(:, nblock) = slater_D(1:nelec)
    R_string(:, :, 1) = slater_R(:, :)
    R_string(:, :, nblock) = slater_R(:, :)

    if(second_trotter) then
      !Q_string(:, 1:nelec, 1) =  matmul(expK_half , Q_string(:, 1:nelec, 1))
      !Q_string(:, 1:nelec, nblock) =  matmul(expK_inv_half,Q_string(:, 1:nelec, nblock))
    end if

    start_block = start_time/ngroup + 1 ! e.g. start_time = 48, ngroup = 5, start_block = 48/5 + 1 = 10
    start_time_by_block = (start_block-1)*ngroup + 1 ! e.g. start_block = 10, ngroup = 5, start_time_by_block = (10-1)*5 + 1 = 46
    !print*,'start_time,start_block,start_time_by_block:',start_time,start_block,start_time_by_block
    if (start_time_by_block > start_time) then
      print *, 'Error in init_g_T0: start_time_by_block should not be larger than start_time'
      stop
    end if

    rmat = Q_string(:, 1:nelec, start_block)

    DO p = start_time_by_block, ntime
      do i_pf = 1, n_phonon_field
        ! left evolve the rmat matrix : expK * expV * Q
        call left_evolve(p, rmat, nelec,i_pf,two_way=.false.)
      end do

      if (mod(p, ngroup) /= 0) cycle; 
      block_num = (p - 1)/ngroup + 2

      call qdr(ns, nelec, rmat, R_string(1:nelec, 1:nelec,block_num), D_string(1:nelec, block_num))

      Q_string(:, 1:nelec, block_num) = rmat
      !print*,'D_string:',block_num,':',D_string(1:nelec,flv,block_num)

    end do

    lmat = conjg(transpose(Q_string(:, 1:nelec, nblock)))
    call cal_g_and_det_T0(lmat,rmat)

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
    complex(dp), allocatable :: lmat(:, :), rmat(:, :), PBP(:, :)
    complex(dp) :: gfast(ns, ns)
    complex(dp) :: ln_cw_fast
    ln_cw_fast = ln_cw
    ln_cw = 0d0
    !print*,'qmat:',qmat

    allocate (rmat(Ns, nelec), lmat(nelec, Ns), PBP(nelec, nelec))

    gfast = g_h(:, :)

    if (forward) then

      block_num = time/ngroup + 1
      rmat = Q_string(:, 1:nelec, block_num - 1)

      DO p = time - ngroup, time - 1
        do i_pf = 1 , n_phonon_field
          call left_evolve(p, rmat, nelec,i_pf,two_way=.false.) !expK * expV * Q
        end do
      end do
      !forall(i = 1:nelec) rmat(:,i) = rmat(:,i) * D_string(i,flv,block_num-1) !Q' = Q*D
      CALL qdr(ns, nelec, rmat, R_string(1:nelec, 1:nelec,block_num), D_string(1:nelec, block_num))
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

      CALL ldq(nelec, Ns, lmat, R_string(1:nelec, 1:nelec,block_num), D_string(1:nelec, block_num))

      Q_string(:, 1:nelec, block_num) = conjg(transpose(lmat))
      D_string(1:nelec, block_num) = D_string(1:nelec, block_num)
      !R_string(1:nelec,1:nelec,flv,block_num) = matmul(R_string(1:nelec,1:nelec,flv,block_num+1),&
      !& R_string(1:nelec,1:nelec,flv,block_num)) ! attention that R_string was not hermitian_conjugated
      rmat = Q_string(:, 1:nelec, block_num - 1)
      !ln_cw = ln_cw + sum(log(D_string(1:nelec,flv,block_num-1)))
    end if
    ! calculate g =  1 - rmat * (lmat*rmat)^(-1) * lmat
    call cal_g_and_det_T0(lmat,rmat)
    ! calculate difference between fast updating g and scratch
    ! difference
    err_fast = maxval(abs(gfast - g_h(:, :)))
    if (err_fast > err_fast_max) then
      print *, 'time:', time, 'error_fast: ', err_fast, 'loop:', loop, 'myid:', myid
      err_fast_max = err_fast
      print *, 'gfast:', maxval(abs(gfast))
      print *, 'g_h:', maxval(abs(g_h(:, :)))
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
    complex(dp), intent(inout) :: mat(Ns, Ns)
    complex(dp) :: gfast(Ns, Ns)
    complex(dp) :: PBP(nelec, nelec)
    complex(dp), intent(inout) :: rmat(Ns, nelec)
    complex(dp) :: lmat(nelec, Ns), D_string_temp(nelec),R_temp(nelec, nelec)
    integer :: i, p, block_num,i_pf
    block_num = time/ngroup + 1 ! e.g. block_num = 46/5 + 1 = 10

    gfast = mat
    DO p = time - ngroup, time - 1
      do i_pf =  n_phonon_field,1,-1
        call left_evolve(p, rmat, nelec,i_pf,two_way=.false.) !expK *  expV * Q
      end do
    end do
    CALL qdr(ns, nelec, rmat, R_temp, D_string_temp)
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

  SUBROUTINE left_evolve(time, mat, n,pf_id,two_way)
    !mat -> expK1*expK2...expKn*mat*exp-Kn...*exp-K2*exp-K1
    IMPLICIT none
    integer, intent(in) :: time, n
    logical, intent(in):: two_way
    complex(dp),intent(inout) :: mat(ns, n)
    integer :: i_pla,pfdim,i,j,n_pla
    integer, intent(in) :: pf_id
    type(pf_data_type), pointer :: p_data
    complex(dp),allocatable :: Q_temp_B(:,:),Q_temp_C(:,:)
    complex(dp),allocatable :: temp_expKV(:,:)
    integer, allocatable :: site_list(:)
    !p_data => ppf%p_data
    !pfdim = ppf%dim
    pfdim = ppf_list(pf_id)%dim
    n_pla = ppf_list(pf_id)%n_plaquette
    allocate(Q_temp_B(pfdim,n),Q_temp_C(pfdim,n), temp_expKV(pfdim,pfdim), site_list(pfdim))
    ! site_list = p_data%pla_site_list

    do i_pla = 1, n_pla
      site_list = pla_site_list(:,i_pla,pf_id) !> get the site list for current phonon field
      Q_temp_B = mat(site_list, :)
      temp_expKV = expkv(:,:,i_pla,pf_id,time)

      !print*,"i_pf,i_pla:",pf_id,i_pla
      !print*,"temp_expKV:",temp_expKV
      !print*,''

      call gemm(Q_temp_C, temp_expKV, Q_temp_B, pfdim, pfdim, n , (1d0,0d0), (0d0,0d0))
      mat(site_list, :) = Q_temp_C
      !!! mat(p_data%pla_site_list(:,i_pla), :) = matmul(p_data%expKV(:,:,i_pla,time) , mat(p_data%pla_site_list(:,i_pla), :))
    end do
    deallocate(Q_temp_B, Q_temp_C, temp_expKV)

    if (two_way) then
      if (n /= ns) stop 'error in left_evolve, dimension not for evolve two-way'
      allocate(Q_temp_B(n,pfdim),Q_temp_C(n,pfdim), temp_expKV(pfdim,pfdim))
        !!! mat(:, p_data%pla_site_list(:,i_pla)) = matmul(mat(:, p_data%pla_site_list(:,i_pla)), p_data%expKV_inv(:,:,i_pla,time))
      do i_pla = 1, n_pla
        site_list = pla_site_list(:,i_pla,pf_id)
        Q_temp_B = mat(:, site_list)
        temp_expKV = expkv_inv(:,:,i_pla,pf_id,time)
        call gemm(Q_temp_C, Q_temp_B, temp_expKV, n, pfdim, pfdim , (1d0,0d0), (0d0,0d0))
        mat(:, site_list) = Q_temp_C
      end do
      deallocate(Q_temp_B, Q_temp_C, temp_expKV)
    end if

    !deallocate(Q_temp_B, Q_temp_C)
  end SUBROUTINE left_evolve

  SUBROUTINE right_evolve(time, mat, n, pf_id,two_way)
    IMPLICIT none
    integer, intent(in) :: time, n
    logical, intent(in) :: two_way
    complex(dp) :: mat(n,ns)
    integer :: i_pla,pfdim,i,j,n_pla
    integer, intent(in) :: pf_id
    type(pf_data_type), pointer :: p_data
    complex(dp),allocatable :: Q_temp_B(:,:),Q_temp_C(:,:)
    complex(dp),allocatable :: temp_expKV(:,:)
    integer, allocatable :: site_list(:) ! (ppf%dim, ppf%n_plaquette)
    pfdim = ppf_list(pf_id)%dim
    n_pla = ppf_list(pf_id)%n_plaquette
    allocate(Q_temp_B(n,pfdim),Q_temp_C(n,pfdim), temp_expKV(pfdim,pfdim), site_list(pfdim))

    do i_pla = 1, n_pla
      site_list = pla_site_list(:,i_pla,pf_id)
      Q_temp_B = mat(:,site_list)
      temp_expKV = expkv(:,:,i_pla,pf_id,time)
      call gemm(Q_temp_C, Q_temp_B, temp_expKV, n, pfdim, pfdim , (1d0,0d0), (0d0,0d0))
      mat(:,site_list) = Q_temp_C
      !!! mat(:,p_data%pla_site_list(:,i_pla)) = matmul(mat(:,p_data%pla_site_list(:,i_pla)),p_data%expKV(:,:,i_pla,time))
    end do
    deallocate(Q_temp_B, Q_temp_C, temp_expKV)

    if (two_way) then

      if (n /= ns) stop 'error in right_evolve, dimension not for evolve two-way'
      allocate(Q_temp_B(pfdim,n),Q_temp_C(pfdim,n), temp_expKV(pfdim,pfdim))
        !!!mat(p_data%pla_site_list(:,i_pla), :) = matmul(p_data%expKV_inv(:,:,i_pla,time),mat(p_data%pla_site_list(:,i_pla), :))
      do i_pla = 1, n_pla
        site_list = pla_site_list(:,i_pla,pf_id)
        Q_temp_B = mat(site_list, :)
        temp_expKV = expkv_inv(:,:,i_pla,pf_id,time)
        call gemm(Q_temp_C, temp_expKV, Q_temp_B,  pfdim, pfdim , n,(1d0,0d0), (0d0,0d0))
        mat(site_list, :) = Q_temp_C
      end do
      deallocate(Q_temp_B, Q_temp_C, temp_expKV)

    end if

  end SUBROUTINE right_evolve

  subroutine cal_g_and_det_T0(lmat,rmat)
    implicit none
    integer :: i
    complex(dp), intent(in) :: lmat(nelec,Ns), rmat(Ns,nelec)
    complex(dp),allocatable :: PBP(:,:)
    complex(dp) :: lmat_T(2*nelec,Ns), rmat_T(Ns,2*nelec)
    complex(dp) :: inner_prod(nelec,nelec)
    complex(dp) :: g_h_buffer(Ns,nelec)
    ln_cw = 0d0
    g_h_buffer = (0d0,0d0)
    g_h = (0d0,0d0)

    allocate(PBP(nelec,nelec))
      !!!PBP = matmul(lmat,rmat)
    call gemm(PBP, lmat, rmat, nelec,  Ns, nelec, (1d0,0d0), (0d0,0d0))
    ln_cw = ln_cw + log((det(nelec, PBP))) !> abs for real only
    ln_cw = ln_cw + sum(log((D_string)))
    call inverse(nelec, PBP) ! PBP^(-1) -> PBP
      !!! g_h = MATMUL(rmat, matmul(PBP, lmat))
    call gemm(g_h_buffer, rmat, PBP, Ns, nelec, nelec, (1d0,0d0), (0d0,0d0))
    call gemm(g_h, g_h_buffer, lmat, Ns, nelec,Ns, (1d0,0d0), (0d0,0d0))

    if(TR_weight) then
      ln_cw = ln_cw + conjg(ln_cw)
    end if
    ln_cw = ln_cw*ncopy
    ln_cw = cmplx(real(ln_cw),mod(aimag(ln_cw),2*pi),kind=dp) ! make sure weight is in the principal branch
  end subroutine cal_g_and_det_T0

  subroutine cal_new_weight_T0()
    ! PQMC version
    ! calculate the relative weight caused by det(), the expKV are  available in pf_data
    ! this subroutine is called in update_global(working version)
    ! the Q,D,R are not changed in this subroutine
    implicit none
    integer flv, p, d, flag, i, i_pf,i_pla
    complex(dp) :: lndet
    complex(dp), allocatable :: tmat(:, :), lmat(:, :), rmat(:, :), expKV_temp(:,:)
    complex(dp), allocatable :: lmat_temp_B(:, :),lmat_temp_C(:,:)
    ln_cw = 0d0
    lndet = 0d0
    d = nelec
    allocate (tmat(d, d), lmat(d, ns), rmat(ns, d))

    lmat = conjg(transpose(Q_string(:, 1:nelec, nblock)))

    do i = 1, d
      lndet = lndet + log((d_string(i, nblock)))
    end do

    flag = 0
    do p = ntime, 1, -1
      do i_pf = n_phonon_field,1,-1
        call right_evolve(p, lmat, nelec,i_pf,two_way=.false.)
      end do
      flag = flag + 1
      if (flag /= ngroup) cycle
      flag = 0
      call zlq(d, ns, lmat, tmat)
      do i = 1, d
        lndet = lndet + log((tmat(i, i)))
      end do
      !print*,'lndet for time',p,'is',lndet
    end do

    rmat = Q_string(:, 1:nelec, 1)

    do i = 1, d
      lndet = lndet + log((d_string(i, 1)))
    end do

    !!!tmat = matmul(lmat, rmat)
    call gemm(tmat, lmat, rmat, d, ns, d, (1d0,0d0), (0d0,0d0))
    lndet = lndet + log((det(d, tmat)))
    !print*,'lndet:',log(abs(det(d,tmat)))
    if(TR_weight) lndet = lndet + conjg(lndet)
    ln_cw = ln_cw + lndet * ncopy
    ln_cw = cmplx(real(ln_cw),mod(aimag(ln_cw),2*pi),kind=dpc) ! make sure weight is in the principal branch
  end subroutine cal_new_weight_T0
! ------------------------------ for finite temperature evolution, not used for now ------------------ --------------

  SUBROUTINE init_g_FT(start_time)

    implicit none
    integer, intent(in) :: start_time
    integer :: start_block, start_time_by_block
    ! recalculate B_string(s) from time = start_time to ntime, yielding green function at time = ntime+1
    ! start_time ranges from 1 to ntime, e.g. start_time = 1 means recalulate all B_strings
    ! start_time = ntime means no evolution, just calculate g_h at ntime+1
    ! start_time and start_block are where the evolution starts.
    integer ::  p, block_num, i_pf,i
    complex(dp) :: R_temp(Ns,Ns)
    complex(dp), allocatable :: lmat(:, :), rmat(:, :)
    IF (.not. allocated(g)) allocate (g(Ns, Ns))
    IF (.not. allocated(g_h)) allocate (g_h(Ns, Ns))
    IF (.not. allocated(Q_string)) allocate (Q_string(Ns, Ns, nblock),D_string(Ns, nblock),R_string(Ns, Ns, nblock))

    ln_cw = 0d0 ! initialize configuration weight
    allocate (rmat(Ns, Ns), lmat(Ns, Ns))

    !! initialize QDR_string(block = 1 and nblock)
    Q_string = 0
    D_string = 1d0
    R_string = 0
    do i = 1 , Ns
      Q_string(i, i, 1) = 1d0
      Q_string(i, i, nblock) = 1d0
      R_string(i, i, 1) = 1d0
      R_string(i, i, nblock) = 1d0
    end do
    if(second_trotter) then
      !Q_string(:, 1:nelec, 1) =  matmul(expK_half , Q_string(:, 1:nelec, 1))
      !Q_string(:, 1:nelec, nblock) =  matmul(expK_inv_half,Q_string(:, 1:nelec, nblock))
    end if

    start_block = start_time/ngroup + 1 ! e.g. start_time = 48, ngroup = 5, start_block = 48/5 + 1 = 10
    start_time_by_block = (start_block-1)*ngroup + 1 ! e.g. start_block = 10, ngroup = 5, start_time_by_block = (10-1)*5 + 1 = 46
    !print*,'start_time,start_block,start_time_by_block:',start_time,start_block,start_time_by_block
    if (start_time_by_block > start_time) then
      print *, 'Error in init_g_T0: start_time_by_block should not be larger than start_time'
      stop
    end if

    rmat = Q_string(:,:, start_block)

    DO p = start_time_by_block, ntime
      do i_pf = 1, n_phonon_field
        ! left evolve the rmat matrix : expK * expV * Q
        call left_evolve(p, rmat, ns,i_pf,two_way=.false.)

      end do

      if (mod(p, ngroup) /= 0) cycle; 
      block_num = (p - 1)/ngroup + 2
      !!! FTQMC is different from here,should not directly call qdr, but need to do QD first then qdr
      do i = 1, Ns
        rmat(:, i) = rmat(:, i) * D_string(i, block_num-1) !Q' = Q*D
      end do

      !call qdr(ns, ns, rmat, R_temp, D_string(:, block_num))
      call udv(ns, ns, rmat, R_temp, D_string(:, block_num))
      ! FTQMC also needs to track the R_string, R_string(:,:,block_num) = matmul(R_string(:,:,block_num), R_string(:,:,block_num-1))
      ! attention that R_string was not hermitian_conjugated
      ! use gemm to do the matrix multiplication for better performance
      call gemm(R_string(:,:,block_num), R_temp, R_string(:,:,block_num-1), ns, ns, ns, (1d0,0d0), (0d0,0d0))

      !!! FTQMC 's difference stops here

      Q_string(:, :, block_num) = rmat
      !print*,'D_string:',block_num,':',D_string(1:nelec,flv,block_num)

    end do
    !!! FTQMC does not need a congj & transpose before Q_string
    !lmat = Q_string(:, :, nblock)
    lmat = conjg(transpose(rmat))

    !call cal_g_and_det_FT(lmat,rmat, D_string(:, nblock), D_string(:, nblock-1), &
    !& R_string(:,:, nblock), R_string(:,:, nblock-1), .true.)
    call cal_g_and_det_FT(lmat,rmat, D_string(:, nblock-1), D_string(:, nblock-1), &
    & conjg(transpose(R_string(:,:, nblock-1))), R_string(:,:, nblock-1), .true.)
    !print*,'max g:',maxval(abs(g_h))

  end SUBROUTINE

 SUBROUTINE cal_g_and_det_FT(Q_l, Q_r, D_l, D_r, R_l, R_r, cal_g)
    implicit none
    ! for FTQMC, g = (1 + B(\tau,0) * B(\beta,\tau))^(-1)
    ! = (1 + Q_r*D_r*R_r * R_l*D_l*Q_l)^(-1)
    ! = (Q_l)^(-1) * [(Q_l Q_r)^(-1) + D_r*(R_r*R_l)*D_l]^(-1) * (Q_r)^(-1)
    ! = (R Q_l)^(-1) * (D)^(-1) * (Q_r Q )^(-1)
    ! where QDR = (Q_l Q_r)^(-1) + D_r * (R_r*R_l) * D_l
    ! ln_cw = log(det(1 +  B(\beta,0))) = log(det(1 + B(\tau,0) * B(\beta,\tau))) = log(det(g^(-1)))
    ! = log(det(Q_r Q)) + log(det(R Q_l)) + sum(log(D)) for FTQMC
    ! use gemm for matrix multiplication for better performance
    integer :: i,j
    logical,intent(in) :: cal_g
    complex(dp), intent(in) :: Q_l(Ns, Ns), Q_r(Ns, Ns), D_l(Ns), D_r(Ns), R_l(Ns, Ns), R_r(Ns, Ns)
    complex(dp) :: Ql_temp(Ns, Ns), Qr_temp(Ns, Ns)
    complex(dp) :: Q_l_Q_r(Ns,Ns), R_r_R_l(Ns,Ns),QDR_temp(Ns,Ns), R_Ql(Ns,Ns), Qr_Q(Ns,Ns)
    complex(dp) :: R_temp(Ns,Ns),D_temp(Ns) ! for QDR calculation
    complex(dp) :: g_buffer(Ns, Ns)
    complex(dp) :: Q_temp(Ns, Ns)
     ! --- Scale Arrays ---
    complex(dp) :: Dl_gt(Ns), Dl_lt(Ns)
    complex(dp) :: Dr_gt(Ns), Dr_lt(Ns)
    ln_cw = 0d0
    ! initialization
    QDR_temp = 0d0
    Q_l_Q_r = 0d0
    R_r_R_l = 0d0
    R_Ql = 0d0
    Qr_Q = 0d0
    g_buffer = 0d0
    Dl_gt = 0d0
    Dl_lt = 0d0
    Dr_gt = 0d0
    Dr_lt = 0d0
    Ql_temp = Q_l
    Qr_temp = Q_r
    ! =========================================================================
    ! STAGE 1: Split the Scales
    ! =========================================================================
    do i = 1, Ns
      ! Split D_l
      if (abs(D_l(i)) > 1.0d0) then
        Dl_gt(i) = D_l(i)
        Dl_lt(i) = (1.0d0, 0.0d0)
      else
        Dl_gt(i) = (1.0d0, 0.0d0)
        Dl_lt(i) = D_l(i)
      end if
      
      ! Split D_r
      if (abs(D_r(i)) > 1.0d0) then
        Dr_gt(i) = D_r(i)
        Dr_lt(i) = (1.0d0, 0.0d0)
      else
        Dr_gt(i) = (1.0d0, 0.0d0)
        Dr_lt(i) = D_r(i)
      end if
    end do



    ! Q_l_Q_r = (Q_l Q_r)^(-1), use that Q_l and Q_r are unitary to avoid matrix inverse
    Qr_temp = conjg(transpose(Q_r))
    Q_l_Q_r = 0d0
    Ql_temp = conjg(transpose(Q_l))
    call gemm(Q_l_Q_r, Qr_temp, Ql_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
    
    ! QDR = R_r*R_l
    call gemm(R_r_R_l, R_r, R_l, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

    ! Scale R_r*R_l by Dr_lt and Dl_lt, and scale (Q_l Q_r)^(-1) by Dl_gt and Dr_gt
    do j = 1, Ns
      do i = 1, Ns
        R_r_R_l(i, j) = R_r_R_l(i, j) * Dr_lt(i) * Dl_lt(j)
        Q_l_Q_r(i, j) = Q_l_Q_r(i, j) / Dr_gt(i) / Dl_gt(j)
      end do
    end do


    ! QDR = (Dr_gt)^(-1) * (Q_l Q_r)^(-1) * (Dl_gt)^(-1) + Dr_lt*(R_r*R_l)*Dl_lt
    QDR_temp = R_r_R_l + Q_l_Q_r

    ! do qdr decomposition to get Q, D and R for QDR
    !call qdr(Ns, Ns, QDR_temp, R_temp, D_temp)
    call udv(Ns, Ns, QDR_temp, R_temp, D_temp)
    ! calculate ln_cw = log(det(Q_r Q)) + log(det(R Q_l)) + sum(log(D))
    Qr_temp = Q_r
    Ql_temp = Q_l
    ln_cw = ln_cw + log(det(Ns, Ql_temp)) + log(det(Ns, Qr_temp)) + sum(log(D_temp)) + sum(log(Dl_gt)) + sum(log(Dr_gt)) &
    & + log(det(Ns,QDR_temp)) + log(det(Ns, R_temp))
    
    ! scale Q_l and Q_r for later use in calculating g
    ! Q_l = Dl_gt * Q_l, Q_r = Q_r * Dr_gt
    !do j = 1, Ns
    !  do i = 1, Ns
    !    Ql_temp(i, j) = Dl_gt(i) * Q_l(i, j)
    !    Qr_temp(i, j) = Q_r(i, j) * Dr_gt(j)
    !  end do
    !end do

    if(cal_g) then
    ! cal g = Ql^(-1) * Dl_gt^(-1) * (R)^(-1) * (D)^(-1) * (Q)^(-1) * Dr_gt^(-1) * Qr^(-1)
     !call inverse(Ns, R_temp)
      R_temp = conjg(transpose(R_temp)) ! use that R is unitary to avoid matrix inverse
      do j = 1, Ns
        do i = 1, Ns
          Ql_temp(i, j) = conjg(Q_l(j, i))
          Qr_temp(i, j) = conjg(Q_r(j, i))
          Q_temp(i, j) = conjg(QDR_temp(j, i)) ! 将 Q_M 的逆存入 Q_temp 以防覆盖
        end do
      end do
    
      ! Ql_temp^(-1) * Dl_gt^(-1)
      do j = 1, Ns
        do i = 1, Ns
          Ql_temp(i, j) = Ql_temp(i, j) / Dl_gt(j)
        end do
      end do
      ! R_Ql = Ql^(-1) * Dl_gt^(-1) * (R)^(-1)
      call gemm(R_Ql, Ql_temp, R_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

      ! R_Ql = R_Ql * (D)^(-1) 
      do j = 1, Ns
        do i = 1, Ns
          R_Ql(i, j) = R_Ql(i, j) / D_temp(j)
        end do
      end do

      !Qr_temp = Dr_gt^(-1) * Qr^(-1)
      do j = 1, Ns
        do i = 1, Ns
          Qr_temp(i, j) = Qr_temp(i, j) / Dr_gt(i)
        end do
      end do
      ! Qr_Q = (Q)^(-1) * Dr_gt^(-1) * Qr^(-1) = (Q)^(-1) * Qr_temp
      call gemm(Qr_Q, Q_temp, Qr_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

      ! g = R_Ql * Qr_Q
      call gemm(g_buffer, R_Ql, Qr_Q, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))
      
      ! 空穴格林函数转换
      g_h(:, :) = - g_buffer(:, :)
      do i = 1, Ns
        g_h(i, i) = g_h(i, i) + 1d0
      end do
    end if

    
    ! make sure weight is in the principal branch
    if(TR_weight) then
      ln_cw = ln_cw + conjg(ln_cw)
    end if
    ln_cw = ln_cw*ncopy
    ln_cw = cmplx(real(ln_cw),mod(aimag(ln_cw),2*pi),kind=dp)

end subroutine


  SUBROUTINE get_g_scratch_FT(time, forward, loop)
    implicit none
    ! for FTQMC
    ! come here when mod(time_count,ngroup)) == 0 including time = ntime
    ! when forward, should give g(time+1)
    ! when backward, should give g(time)
    integer, intent(in) :: time, loop
    logical, intent(in) :: forward
    integer :: i,j, p, block_num,i_pf
    complex(dp) :: lmat(Ns,Ns), rmat(Ns,Ns), R_temp(Ns,Ns), D_temp(Ns)
    complex(dp) :: gfast(Ns, Ns)
    complex(dp) :: ln_cw_fast
    ln_cw_fast = ln_cw
    ln_cw = 0d0
    R_temp = 0d0
    D_temp = 0d0
    !print*,'qmat:',qmat
    gfast = g_h(:, :)

    if (forward) then

      block_num = time/ngroup + 1
      rmat = Q_string(:, :, block_num - 1)

      DO p = time - ngroup, time - 1
        do i_pf = 1 , n_phonon_field
          call left_evolve(p, rmat, Ns,i_pf,two_way=.false.) !expK * expV * Q
        end do
      end do
      !Q' = Q*D
      do i = 1, Ns
        rmat(:, i) = rmat(:, i) * D_string(i, block_num-1)
      end do
      CALL udv(Ns, Ns, rmat, R_temp, D_string(:, block_num))
      Q_string(:, :, block_num) = rmat
      call gemm(R_string(:,:,block_num), R_temp, R_string(:,:,block_num-1), Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

      lmat = Q_string(:, :, block_num + 1)
      !update ln_cw

      call cal_g_and_det_FT(lmat,rmat, D_string(:, block_num+1), &
    & D_string(:, block_num), R_string(:,:, block_num+1), R_string(:,:, block_num),.true.)
    else !evolve backwards

      block_num = time/ngroup + 2  ! depends on size of ngroup
      !print*,"time,block_num:",time,block_num
      if (ngroup == 1) block_num = block_num - 1
      lmat = Q_string(:, :, block_num + 1)
      !print*,"lmat:",real(lmat)
      DO p = time + ngroup - 1, time, -1
        do i_pf = n_phonon_field,1,-1
          call right_evolve(p, lmat, Ns,i_pf,two_way=.false.)!Q*expK * expV
        end do
      end do
      !print*,"lmat:",real(lmat)
      ! forall(i = 1:nelec) lmat(i,:) = lmat(i,:) * D_string(i,flv,block_num+1) !Q' = Q*D
      do i = 1, Ns
        lmat(i, :) = lmat(i, :) * D_string(i, block_num+1)
      end do
      !print*,"lmat:",real(lmat)
      CALL udv(Ns, Ns, lmat, R_temp, D_string(:, block_num))
      !print*,"lmat:",real(lmat)
      Q_string(:, :, block_num) = lmat
      call gemm(R_string(:,:,block_num), R_string(:,:,block_num+1),R_temp, Ns, Ns, Ns, (1d0,0d0), (0d0,0d0))

      rmat = Q_string(:, :, block_num - 1)
      !ln_cw = ln_cw + sum(log(D_string(1:nelec,flv,block_num-1)))

      ! calculate g =  (1 + B(\tau,0) * B(\beta, \tau))^(-1)
      call cal_g_and_det_FT(lmat,rmat, D_string(:, block_num), &
      & D_string(:, block_num-1), R_string(:,:, block_num), R_string(:,:, block_num-1),.true.)
    end if

    ! calculate difference between fast updating g and scratch
    ! difference
    err_fast = maxval(abs(gfast - g_h(:, :)))
    if (err_fast > err_fast_max) then
      print *, 'time:', time, 'error_fast: ', err_fast, 'loop:', loop, 'myid:', myid
      err_fast_max = err_fast
      print *, 'gfast:', maxval(abs(gfast))
      print *, 'g_h:', maxval(abs(g_h(:, :)))
      !stop 'Error in update: difference too large'
    end if
    err_fast = abs((real(ln_cw - ln_cw_fast)))
    if (err_fast > err_fast_max) then
      err_fast_max = err_fast
      print *, 'time:', time, 'error_fast_det: ', err_fast, 'loop:', loop, 'myid:', myid
      print *, 'lnw_fast,lnw_from_scratch:', ln_cw_fast, ln_cw
      !stop 'Error in update lnw: difference too large'
    end if
  end SUBROUTINE get_g_scratch_FT

  subroutine cal_new_weight_FT()
    ! FTQMC version
    ! calculate the relative weight caused by det()
    ! this subroutine is called in update_global(working version)
    ! the Q_string,D_string,R_string are intact in this subroutine
    implicit none
    integer flv, p, flag, i, i_pf,i_pla,block_num
    complex(dp) :: lndet
    complex(dp), allocatable :: lmat(:, :), rmat(:, :),d_temp(:),R_temp(:,:),R_new(:,:),R_temp2(:,:)
    ln_cw = 0d0
    lndet = 0d0
    allocate (lmat(Ns, ns), rmat(ns, Ns),d_temp(Ns),R_temp(Ns,Ns),R_temp2(Ns,Ns),R_new(Ns,Ns))
    rmat = Q_string(:,:, 1)
    d_temp = D_string(:,1)
    R_new = R_string(:,:,1)

    DO p = 1, ntime
      do i_pf = 1, n_phonon_field
        ! left evolve the rmat matrix : expK * expV * Q
        call left_evolve(p, rmat, ns,i_pf,two_way=.false.)

      end do

      if (mod(p, ngroup) /= 0) cycle; 
      block_num = (p - 1)/ngroup + 2
      !!! FTQMC is different from here,should not directly call qdr, but need to do QD first then qdr
      do i = 1, Ns
        rmat(:, i) = rmat(:, i) * d_temp !Q' = Q*D
      end do

      call udv(ns, ns, rmat, R_temp, d_temp)
      ! FTQMC also needs to track the R_string, R_string(:,:,block_num) = matmul(R_string(:,:,block_num), R_string(:,:,block_num-1))
      ! attention that R_string was not hermitian_conjugated
      ! use gemm to do the matrix multiplication for better performance
      R_temp2 = R_new
      call gemm(R_new, R_temp, R_temp2, ns, ns, ns, (1d0,0d0), (0d0,0d0))

    end do
    !!! FTQMC does not need a congj & transpose before Q_string
    lmat = Q_string(:, :, nblock)
    call cal_g_and_det_FT(lmat,rmat, D_string(:, nblock), d_temp, R_string(:,:, nblock), R_new, .false.)

  end subroutine cal_new_weight_FT

end module evolution