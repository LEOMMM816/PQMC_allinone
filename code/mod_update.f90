module update

  use ran_state
  use evolution
  implicit none

contains

!======================================================================!
!                             LOCAL UPDATE                             !
!======================================================================!

  subroutine update_phonon_field(ppf, time)

    !> update phonon field, calculate exp(delta E), if accepted, update Kmat, green function
    !> this subroutine is called in Main_PQMC.f90
    !> ppf is the phonon field type, time is the time step
    !> note that boson_field and the expKV and expKV_inv have been updated
    !> update ph_field, calculate exp(delta E), if accepted, update Kmat, green function
    implicit none
    integer, intent(in) :: time
    type(pf_type), intent(in) :: ppf
    type(pf_data_type), pointer :: p_data
    complex(dp) :: R,R_elec ! relative boltzmann weight &
    real(dp) ::  delta_energy, bf_jump, delta_energy_elec
    real(dp) :: old_bf,prob
    complex(dp), dimension(ppf%dim,ppf%dim) :: expKV_old,expK_inv_old
    complex(dp),dimension(ppf%dim,ppf%dim) :: Delta_sub,  R_mat_sub, R_temp_sub, g_h_sub,g_sub
    complex(dp), dimension(Ns, Ns) :: g_temp
    complex(dp) :: gh_part_nd(Ns,ppf%dim) !> Ns*d of G_h,auxillary matrix to construct updated g
    complex(dp) :: g_part_dn(ppf%dim,Ns)!> d*Ns of G,auxillary matrix to construct updated g
    complex(dp) :: gh_part_nd_temp(Ns,ppf%dim),scalar_val
    integer ::  i,j,k,i_pla
    integer :: site_list(ppf%dim)
    logical :: flag = .true., update_accept
    real(dp), pointer :: p_bf
    integer :: pfdim
    !if(.not.allocated(g_debug)) allocate(g_debug(Ns,Ns))

    !initialization
    p_data => ppf%p_data
    Delta_sub = 0d0
    g_h_sub = 0d0
    R_mat_sub = 0d0
    g_temp = 0d0
    g_sub = 0d0
    !g_debug = 0d0
    gh_part_nd = 0d0
    gh_part_nd_temp = 0d0
    g_part_dn = 0d0
    delta_energy = 0d0
    delta_energy_elec = 0d0
    bf_jump = 0d0
    R_temp_sub = 0d0
    pfdim = ppf%dim !> dimension of phonon field
    flag = .true.

      !! fix two endpoints
    ! if ( time ==1 .or. time == ntime ) return

      !! preparation

    do i_pla = 1 , ppf%n_plaquette

      ! get the site index of the plaquette
      site_list = p_data%pla_site_list(:,i_pla)
      ! get the old field value
      p_bf => boson_field(p_data%bf_list(i_pla), time)
      old_bf = p_bf
      ! get g_sub and g_h_sub
      g_h_sub = g_h(site_list, site_list)
      g_sub = -g_h_sub
      do j = 1, pfdim
        g_sub(j, j) = g_sub(j, j) + 1d0
      end do
      update_accept = .false.
      do i = 1, n_local_update
        R = 1d0
        R_elec = 1d0

        ! new bf_jump value
        bf_jump = jump_distance*ran_sysm()

        ! calculate energy difference due to boson field: delta_energy
        call cal_energy_difference(p_data%bf_list(i_pla), time, bf_jump, delta_energy)
        delta_energy = -delt*delta_energy
        R = R * exp(delta_energy)

        ! update boson field
        p_bf = p_bf + bf_jump

        ! calculate det(B')/det(B), energy difference due to electrons
        Delta_sub = calculate_Delta(ppf,i_pla,time,p_bf)
        !> R_sub = I  + G_h_sub * Delta_sub
        R_mat_sub = matmul(g_h_sub,Delta_sub)
        do j = 1 , pfdim
          R_mat_sub(j, j) = R_mat_sub(j, j) + 1d0
        end do
        R_elec = R_elec * det(pfdim,R_mat_sub)! R_elec = det(I + delta*G_h)
        
        IF(TR_weight) R_elec = R_elec * CONJG(R_elec) ! for TR system, the weight is squared
        R_elec = R_elec**ncopy

        ! total weight ratio
        R = R*(R_elec)
        !print*, 'i_pla,time,bf_jump,delta_energy,R_elec,R:',i_pla,time,bf_jump,delta_energy,R_elec,R
        ! calculate the acceptance probability
        call cal_prob(log(R), prob,"local update")

        ! acceptance test
        !if(.true.) then
        if ((prob < rands())) then
          ! reject
          p_bf = p_bf - bf_jump
        else
          update_accept = .true.
          ! update pf_datam
          call get_expKV(p_data%expKV(:,:,i_pla,time),ppf,i_pla,p_bf,inv=.false.)
          call get_expKV(p_data%expKV_inv(:,:,i_pla,time),ppf,i_pla,p_bf,inv=.true.)

          ! update G_sub
          ! g_h_sub = g_h_sub + g_h_sub * Delta * R^(-1) * (I - g_h_sub)
          call inverse(pfdim, R_mat_sub) ! R_mat_sub => R_mat_sub^(-1)
          R_temp_sub = matmul(Delta_sub,R_mat_sub)
          ! g_h_sub = g_h_sub + matmul(matmul(g_h_sub,R_temp_sub),g_sub)
          g_sub = matmul(R_temp_sub, g_sub)
          g_sub = matmul(g_h_sub, g_sub)
          g_h_sub = g_h_sub + g_sub
          g_sub = -g_h_sub
          do j = 1 , pfdim
            g_sub = g_sub + 1d0
          end do
          ln_cw = ln_cw + log(real(R_elec))
        end if
      end do  ! end of local update for one plaquette

      if (update_accept) then
        !if(.false.) then
         !! final update ph & g
        positive_accept = positive_accept + 1
        ! update green functiuon, the derivation see my note called "delayed update"
        ! note that boson_field and the expKV and expKV_inv have been updated
        g_h_sub = g_h(site_list, site_list)
        gh_part_nd = g_h(:, site_list)
        do i = 1, pfdim
          g_part_dn(i,:) = - g_h(site_list(i),:)
          g_part_dn(i, site_list(i)) = g_part_dn(i, site_list(i)) + 1d0
        end do
        call get_expKV(Delta_sub,ppf,i_pla, old_bf, inv=.true.)
        ! now Delta_sub = exp(delt*KV), the old one; p_data%expKV contains the new ones: exp(-delt*KV')
        ! we need to do: Delta_sub = expKV * exp(-delt*KV') - I or Delta_sub = Delta_sub * p_data%expKV - I
        Delta_sub  = matmul(Delta_sub, p_data%expKV(:,:,i_pla,time))
        do j = 1 , pfdim
          Delta_sub(j, j) = Delta_sub(j, j) - 1d0
        end do
        ! R_mat_sub = (I + g_h_sub * Delta_sub)
        R_mat_sub = matmul(g_h_sub,Delta_sub)
        do i = 1 , pfdim
          R_mat_sub(i, i) = R_mat_sub(i, i) + 1d0
        end do
        ! R_mat_sub => R_mat_sub^(-1)
        call inverse(pfdim,R_mat_sub)
        ! R_mat_sub  => Delta_sub * R_mat_sub^(-1)
        R_mat_sub = matmul(Delta_sub,R_mat_sub)

        ! update g_h and g, g_h = g_h + g_h * R_mat_sub * g
        gh_part_nd_temp = matmul(gh_part_nd,R_mat_sub)
        !DO j = 1, Ns
        !  DO k = 1, pfdim
        !    scalar_val = g_part_dn(k, j)
        !    ! 最内层循环遍历 i (N)，内存连续，这是向量化的黄金位置
        !    DO i = 1, Ns
        !      g_h(i, j) = g_h(i, j) + gh_part_nd_temp(i, k) * scalar_val
        !    END DO
        !  END DO
        !END DO
        call gemm(g_h, gh_part_nd_temp, g_part_dn, Ns,  pfdim, Ns,(1d0,0d0), (1d0,0d0))

        ! finally, get g from g_h

        
        ! update Kmat
        !call update_Kmat(bond,time)
      else
        p_bf = old_bf
        negative_accept = negative_accept + 1
        ! if not accepted, restore the old boson field
        !print*,'non-accepted ln_cw:',ln_cw
      end if
    end do ! end of loop over plaquettes

  end subroutine

  Subroutine cal_energy_difference(bf_index, time, bf_jump, delta_energy)
    implicit none
    integer ::  time,bf_index
    real(dp) :: bf_jump, delta_energy, E_1, E_2, ph_be, ph_af, ph_now, ph_new
    ph_be = boson_field(bf_index, modi(time - 1, ntime))
    ph_af = boson_field(bf_index, modi(time + 1, ntime))
    ph_now = boson_field(bf_index, time)
    ph_new = ph_now + bf_jump
    E_1 = 0.5*M*((ph_be - ph_now)**2 + (ph_af - ph_now)**2)/(delt)**2 + 0.5*D*ph_now**2 - biased_phonon * ph_now
    E_2 = 0.5*M*((ph_be - ph_new)**2 + (ph_af - ph_new)**2)/(delt)**2 + 0.5*D*ph_new**2 - biased_phonon * ph_new
    delta_energy = (E_2 - E_1)
    if(.false.) then
      print*, 'ph_be,ph_now,ph_af:',ph_be,ph_now,ph_af
      print*, 'Ek_1:',0.5*M*((ph_be - ph_now)**2 + (ph_af - ph_now)**2)/(delt)**2
      print*, 'Ep_1:',0.5*D*ph_now**2
      print*, 'E_1:',E_1
      print*, 'ph_be,ph_new,ph_af:',ph_be,ph_new,ph_af
      print*, 'Ek_2:',0.5*M*((ph_be - ph_new)**2 + (ph_af - ph_new)**2)/(delt)**2
      print*, 'Ep_2:',0.5*D*ph_new**2
      print*, 'E_2:',E_2
      print*, 'delta_energy:',delta_energy

      stop
    end if
  end subroutine

  FUNCTION calculate_Delta(ppf,i_pla, time,bf_new) result(re)
    implicit none
    integer, intent(in)::  i_pla,time! n is dimension of field
    type(pf_type), intent(in) :: ppf
    complex(dp) :: re(ppf%dim, ppf%dim)
    real(dp),pointer ::  bf_new
    real(dp) :: bf_value
    integer :: i,j
    bf_value = bf_new
    call get_expKV(re,ppf,i_pla,bf_value,inv = .false.) ! re = expKV
    re = matmul(ppf%p_data%expKV_inv(:,:,i_pla,time),re) ! delta_sub = exp(delt*KV) * exp(-delt*KV') - I
    do i = 1 , ppf%dim
      re(i, i) = re(i, i) - 1d0
    end do
    ! p_Data%expKV_inv(:,:,i_pla,time) is exp(delt*KV)
  end FUNCTION

!======================================================================!
!                            GLOBAL UPDATE                             !
!======================================================================!

  subroutine update_global()
    implicit NONE
    real(dp) :: oldfield(Ns, ntime)
    complex(dp) :: oldconfigurationweight
    real(dp) :: ph_ln_cw,prob,test_ln_cw
    type(pf_data_type) :: old_pdata(n_phonon_field)
    integer :: i_pf
    oldfield = boson_field
    oldconfigurationweight = ln_cw
    old_pdata = pf_data_list
    ln_cw = 0d0
    test_ln_cw = 0d0
    call generate_newfield_global() ! generate new pf and sum over the weight caused by energy difference
    !call generate_newfield_space(bf_jump, global_update_Kvec, global_update_distance)
    ph_ln_cw = real(ln_cw)

    !call cal_new_weight() ! calculate the relative weight caused by det()
    !test_ln_cw = real(ln_cw)
    call cal_new_weight_gemm() ! calculate the relative weight caused by det()
    !print*,'ln_cw:',ln_cw
    !print*,'ph_ln_cw:',ph_ln_cw
    !print*,'ln_cw + ph_ln_cw - oldconfigurationweight:',ln_cw + ph_ln_cw - oldconfigurationweight

    call cal_prob(ln_cw + ph_ln_cw - oldconfigurationweight,prob,"global update")

    if (ran() < prob) then ! accept this global update
      global_accept = global_accept + 1
      updated = .true.
      ! need to update expKV and expKV_inv in pf_list
      ! do it in generate_newfield_global subroutine
      !print*,'accepted!'
    else
      boson_field = oldfield
      ln_cw = oldconfigurationweight
      pf_data_list = old_pdata
      !Kmat = old_Kmat
      global_reject = global_reject + 1
    end if

  end subroutine

  subroutine generate_newfield_global()
    implicit none
    integer :: i, site, time, i_pla
    integer, dimension(n_global_update) :: updated_ind !
    real(dp) :: bf_temp_ntime(ntime), bf_jump, bf_temp_space(n_boson_field)
    real(dp) :: delta_energy,K_vec(Lat%dim), incell_phase(bf_sets)
    real(dp) :: old_field(n_boson_field, ntime)
    old_field = boson_field
    if((.not.global_update_exchange) .and. (.not. Kspace_GU)) then
      !! new field by site
      do i = 1, n_global_update
        updated_ind(i) = irands(n_boson_field) + 1
        if (global_update_flip) then
          bf_temp_ntime = - boson_field(updated_ind(i), :)
        elseif (global_update_shift) then
          bf_jump = ran_sysm()* global_update_distance
          bf_temp_ntime = boson_field(updated_ind(i), :) + bf_jump
        end if
        !bf_temp_ntime(1) = end_field
        !bf_temp_ntime(ntime) = end_field
      !! energy difference
        delta_energy = 0d0
        if (.not. global_update_flip) then
          delta_energy = 0.5*D*sum(bf_temp_ntime**2 - boson_field(updated_ind(i), :)**2)
        end if
        delta_energy = delta_energy + sum(- biased_phonon * (bf_temp_ntime - boson_field(updated_ind(i), :)))
        boson_field(updated_ind(i), :) = bf_temp_ntime
        ln_cw = ln_cw - delta_energy*delt
      end do
    elseif(global_update_exchange) then!! exchange nearby ph fields
      ! no need to calculate energy difference since it is zero
      do i = 1, n_global_update
        delta_energy = 0d0
        updated_ind(i) = irands(n_boson_field) + 1
        site = irands(n_boson_field) + 1
        bf_temp_ntime = boson_field(updated_ind(i), :)
        boson_field(updated_ind(i), :) = boson_field(site, :)
        boson_field(site, :) = bf_temp_ntime
      end do
    end if

    if(Kspace_GU) then
      K_vec = 2*PI*[0.5d0,0.5d0]  ! Initialize K_vec with appropriate values
      ! incell_phase = [1d0, -1d0]  ! Initialize incell_phase with appropriate values
      incell_phase(1) = irands(2)*2d0 - 1d0 ! random phase + or - 1
      incell_phase(2) = irands(2)*2d0 - 1d0
      ! generate new field in K space
      call generate_newfield_space(bf_temp_space, K_vec, incell_phase)
      do time = 1, ntime
        do i = 1,n_boson_field
          boson_field(i,time) = boson_field(i,time) + bf_temp_space(i)
        end do
      end do
      !boson_field(:, 1) = end_field
      !boson_field(:, ntime) = end_field
      ! calculate the energy difference
      delta_energy = 0d0
      do i = 1, n_boson_field
        delta_energy = delta_energy + 0.5*D*sum(boson_field(i, :)**2 - old_field(i, :)**2)
        delta_energy = delta_energy + sum(- biased_phonon * (boson_field(i, :) - old_field(i, :)))
      end do
      ln_cw = ln_cw - delta_energy*delt
    end if

    ! update the expKV and expKV_inv in pf_list
    do i = 1, n_phonon_field
      do time = 1, ntime
        do i_pla = 1, ppf_list(i)%n_plaquette

          if(pf_list(i)%V_exist) then
            call get_expKV(pf_list(i)%p_data%expKV(:,:,i_pla,time),ppf_list(i),i_pla, &
            & boson_field(ppf_list(i)%p_data%bf_list(i_pla), time), inv=.false.)
            call get_expKV(pf_list(i)%p_data%expKV_inv(:,:,i_pla,time),ppf_list(i),i_pla, &
            & boson_field(ppf_list(i)%p_data%bf_list(i_pla), time), inv=.true.)
          end if
        end do
      end do
    end do
    !call set_kmat()
  end subroutine

  subroutine generate_newfield_space(bf_jump, K_vec, incell_phase)
    implicit none
    real(dp), intent(in) :: incell_phase(bf_sets)
    real(dp), intent(in) :: K_vec(Lat%dim) !K_vec's components belong to [0,2pi]
    real(dp), intent(out) :: bf_jump(n_boson_field)
    real(dp) :: phi, amplitude,distance
    integer :: i_bf,i_cell
   ! print*,"interphase:", incell_phase
    bf_jump = 0d0
    distance = 2 * global_update_distance/lat%N_cell ! can be adjusted
    amplitude =  distance * rands()
    phi = 2d0*pi*rands()
    do i_cell = 1, Lat%N_cell
      do i_bf = 1, bf_sets
        bf_jump(bf_sets*(i_cell-1) + i_bf) = &
        & incell_phase(i_bf) * amplitude*(cos(sum(K_vec * p_cells(i_cell)%rpos) + phi))
      end do
    end do
  end subroutine

  subroutine cal_new_weight()
    ! calculate the relative weight caused by det(), the expKV are not available
    ! this subroutine is called in update_global(backup version)
    ! note that boson_field have been updated but the expKV and expKV_inv need to be calculated explicitly
    ! the Q,D,R are not changed in this subroutine
    implicit none
    integer flv, p, d, flag, i, i_pf,i_pla
    complex(dp) :: lndet
    complex(dp), allocatable :: tmat(:, :), lmat(:, :), rmat(:, :), expKV_temp(:,:)
    complex(dp),allocatable :: lmat_temp_B(:, :),lmat_temp_C(:,:)
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
        allocate (expKV_temp(ppf_list(i_pf)%dim, ppf_list(i_pf)%dim), & 
        & lmat_temp_B(d,ppf_list(i_pf)%dim), lmat_temp_C(d,ppf_list(i_pf)%dim))
        do i_pla = 1, ppf_list(i_pf)%n_plaquette
          if(pf_list(i_pf)%V_exist) then
            call get_expKV(expKV_temp, ppf_list(i_pf), i_pla, boson_field(ppf_list(i_pf)%p_data%bf_list(i_pla), p), inv=.false.)
          else
            expKV_temp = pf_list(i_pf)%p_data%expKV(:,:,i_pla,p) ! if no V, expKV = I
          end if
          
          ! right evolve the lmat matrix : expK * expV * Q
          lmat(:,ppf_list(i_pf)%p_data%pla_site_list(:,i_pla)) = & 
          & matmul(lmat(:,ppf_list(i_pf)%p_data%pla_site_list(:,i_pla)),expKV_temp)
        end do
        deallocate(expKV_temp, lmat_temp_B, lmat_temp_C)
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

    tmat = matmul(lmat, rmat)
    
    lndet = lndet + log((det(d, tmat)))
    !print*,'lndet:',log(abs(det(d,tmat)))
    if(TR_weight) lndet = lndet * conjg(lndet)
    ln_cw = ln_cw + lndet * ncopy
    ln_cw = complex(real(ln_cw),mod(aimag(ln_cw),2*pi)) ! make sure weight is in the principal branch
  end subroutine

  subroutine cal_new_weight_gemm()
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
        call right_evolve(p, lmat, nelec,ppf_list(i_pf),two_way=.false.)
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
    if(TR_weight) lndet = lndet * conjg(lndet)
    ln_cw = ln_cw + lndet * ncopy
    ln_cw = complex(real(ln_cw),mod(aimag(ln_cw),2*pi)) ! make sure weight is in the principal branch
  end subroutine
!======================================================================!
!                            Metropolis prob                         !
!======================================================================!

  subroutine cal_prob(w, prob,sign_err)
    implicit none
    complex(dp), intent(in) :: w
    real(dp), intent(out) :: prob
    character(len=*) :: sign_err
    real(dp) :: R_temp
    complex(dp) :: weight
    ! weight is the energy differece between two configurations
    ! prob is the acceptance probability
    ! R should be exp(weight)
    weight = w
    weight = complex(real(weight),mod(aimag(weight),2 * pi)) ! make sure weight is in the principal branch
    if(abs(aimag(weight)) > 1d-2 .and. abs(abs(aimag(weight)) - 2*pi) > 1d-2) then
      print *, "SIGN PROBLEM,weight:", weight
      print*, "happens in  ", trim(adjustl(mpi_info))//sign_err
      stop
    end if
    R_temp = exp(min(10d0, real(weight)))
    if (Metropolis) then
      prob = min(1d0, real(R_temp))
      if (greatest_decent) then
        prob = floor(prob)
      end if
    else!heat bath
      prob = real(R_temp/(1 + R_temp)) ! accept probobility
    end if
  end subroutine

end module update