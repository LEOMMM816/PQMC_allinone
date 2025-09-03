module zxy_measure


  implicit none


! measurement parameter
  character*30 :: ci, ci2 ! name of file
  character*100 :: str
  integer :: obs, tau, site_1,site_2
  integer :: head = 1, bin = 1!> the current position in current bin and the current bin
  integer :: bin_tau = 1,head_tau = 1 !> especially for tau
  integer :: raw_output_list(4) = (/4,12,30,58/)
  real(8) :: phase_factor
! matrix declare
  complex(kind=8), allocatable :: observables(:, :, :) ! (n_obs,one_bin,nbin_per_core)
  complex(kind=8), allocatable :: aver_bin(:, :) ! (nbin_per_core,n_obs)
  complex(kind=8), allocatable :: obs_aver(:), obs_err(:) ! occupation number and its variance
  complex(kind=8), allocatable :: corfk_phase(:,:) ! phase factor of the (corf_length,corf_length) (should be complex in general)
  complex(8),allocatable :: corf_bin(:,:,:), corfk_bin(:,:,:)! corf_bin function (length, ,n_corf , nbin_per_core)
  integer,allocatable :: corf_s2_list(:,:), corf_s1(:), corf_s2_next(:)
  complex(8),allocatable :: corf_aver(:,:), corf_err(:,:) ! (corf_length,n_corf)
  complex(8),allocatable :: corfk_aver(:,:), corfk_err(:,:) ! (corf_length,n_corf)
! displaced time green funcitons
  complex(8), allocatable :: g_h_tau(:,:,:),g_tau(:,:,:) ! (Ns,Ns,ntau+1) equal-time gf at all tau
  complex(8), allocatable :: g_h_dtau(:,:,:),g_dtau(:,:,:) ! (Ns,Ns,ntau) displaced-time gf at all (tau,tau+1)
  complex(8), allocatable :: g_h_tt(:,:,:),g_tt(:,:,:) ! (Ns,Ns,ntau+1) displaced-time gf at (0,tau), tau from 0 to ntau
  complex(8), allocatable :: ctau_bin(:,:,:),ctau_aver(:,:),ctau_err(:,:) ! (ntau+1,n_ctau,nbin_per_core)
  complex(8), allocatable :: cw_bin(:,:,:),cw_aver(:,:),cw_err(:,:) ! (ntau+1,n_ctau,nbin_per_core)
  complex(8), allocatable :: ctau_phase(:,:) ! (ntau+1,ntau+1)
! observeable declare
  type :: observable_type
    character*30 :: name
    integer :: obs_op_list(4,2,4)
    complex(8) :: para(4)
  end type

 

CONTAINS

  subroutine init_meas()
    implicit none
    integer :: s1, s2 , i , n
    real(8) :: wave_vector ! phase factor of the position (should be complex in general)
    if (.not. allocated(observables)) then
      allocate (observables(n_obs,one_bin, nbin_per_core))
      allocate (corf_bin(corf_length, n_corf,nbin_per_core),corfk_bin(corf_length,n_corf,nbin_per_core))
      allocate (corf_s2_list(corf_length,4),corfk_phase(corf_length,corf_length))
      allocate (corf_s1(4),corf_s2_next(corf_length))
      allocate (ctau_phase(ntau+1,ntau+1))
      allocate (ctau_bin(ntau+1,n_ctau,nbin_per_core),ctau_aver(ntau+1,n_ctau),ctau_err(ntau+1,n_ctau))
      allocate (cw_bin(ntau+1,n_ctau,nbin_per_core),cw_aver(ntau+1,n_ctau),cw_err(ntau+1,n_ctau))
      La = Lat%dlength(1)
      ctau_bin = 0d0
      ctau_aver = 0d0
      ctau_err = 0d0
      observables = 0d0 !
      corf_bin = 0d0
      corf_s2_list = 0
      corfk_bin = 0d0
      corfk_phase = 0d0
      corf_s1 = 0
      corf_s2_next = 0
      ctau_phase = 0d0
    end if

      !! calculate the phase factor of the fourier transformation for each k, number of k is corf_length
    do s1 = 1, corf_length
      ! s1  is the index of k , - pi to (corf_length - 1)/La * pi

      do s2 = 1, corf_length
        ! s2 is the index of r, 0 to (corf_length - 1)
        corfk_phase(s1, s2) = exp(complex(0d0, 2d0 * pi * (s1 - 1-La/2) * (s2 - 1) / La))
      end do
    end do
      !! calculate the phase factor of the fourier transformation for each w, number of w is ntau+1
    do s1 = 1, ntau + 1
      ! s1  is the index of w , - pi to (ntau + 1)/La * pi
      wave_vector =  2 * pi * (s1-1) / (ntau + 1) - pi + pi/(ntau + 1)
      !wave_vector =  2 * pi * (s1-1) / (ntau) - pi
      do s2 = 1, ntau + 1
        ! s2 is the index of t, 0 to (ntau + 1)
        ctau_phase(s1, s2) = 2 * exp(complex(0d0, wave_vector * (s2-1)))
      end do
      ctau_phase(s1, 1) = ctau_phase(s1, 1)*0.5d0
    end do

        
    if(.not. allocated(obs_name)) then
      allocate (obs_name(n_obs))
      obs_name(1) = 'ph_PE'
      obs_name(2) = 'ph_KE'
      obs_name(3) = 'e_KEx'
      obs_name(4) = 'e_KEy'
      obs_name(5) = 'e_PE'
      obs_name(6) = 'e_ph_E'
      obs_name(7) = 'e_TE'
      obs_name(8) = 'ph_X'
      obs_name(9) = 'N_tot'
      obs_name(10) = 'Dou_occ'
      obs_name(11) = 'acc_global'
    end if
    if(.not. allocated(corf_name)) then
      allocate (corf_name(n_corf))
      n  = n_suit_corf
      corf_name(1) = 'Ge'
      corf_name(2) = 'Gpair'
      corf_name(3) = 'nu-nu'
      corf_name(4) = 'nu-nd'
      corf_name(5) = 'den-den'
      corf_name(6) = 'X-X'
      corf_name(7) = 'VBS'
      corf_name(8) = 'JJ'
      do i =1 , n
        corf_name(1 * n + i) = trim(corf_name(i)) // '_AB'
        corf_name(2 * n + i) = trim(corf_name(i)) // '_BA'
        corf_name(3 * n + i) = trim(corf_name(i)) // '_BB'
        corf_name(i) = trim(corf_name(i)) // '_AA'
      end do
      corf_name(n * 4 + 1) = 'JJT'
    end if

    if(.not. allocated(ctau_name)) then
      allocate (ctau_name(n_ctau))
      ctau_name(1) = 'JJ'
      n = n_suit_ctau
      do i =1 , n
        ctau_name(1 * n + i) = trim(ctau_name(i)) // '_AB'
        ctau_name(2 * n + i) = trim(ctau_name(i)) // '_BA'
        ctau_name(3 * n + i) = trim(ctau_name(i)) // '_BB'
        ctau_name(i) = trim(ctau_name(i)) // '_AA'
      end do
      ctau_name(n * 4 + 1) = 'JJT'
    end if
  end subroutine init_meas

  subroutine take_measurement(time, forward)
    implicit none
    integer, intent(in) :: time
    logical, intent(in) :: forward
    integer ::  c1,c2,s1, s2, s1_next! used in fourier transformation
    complex(kind=8) :: factor ! phase factor of the position (should be complex in general)
    complex(8) :: obs_site,  t_s1_c, t_s2_c
    real(8) :: t_s1
    real(8),allocatable :: t_s2(:)
    integer :: i,j,corf_offset,corf_s1_next
    complex(8),allocatable :: corf_temp(:)
    type(cell_type), pointer :: pc1,pc2
    ! 2nd trotter
    if(second_trotter) then

    end if
    if (bin > nbin_per_core) then
      !    print*,'MCS: ' ,loop
      stop 'error in measurement: number of measurements overwhelms bins.'
    end if

    if(.not.allocated(t_s2)) allocate (t_s2(corf_length))
    if(.not.allocated(corf_temp)) allocate (corf_temp(corf_length))
      !! MEASUREMENT BEGINS

      !! s1,s2 can be treated as (s1,s2) or (s1,s1+s2), later when (s1,s2) plays the role of coordinates

      !! obs1: phonon potential energy
    obs_site = 0d0
    obs_site = obs_site + 0.5d0 * D * sum(boson_field(:,time)**2) - biased_phonon * sum(boson_field(:,time))
    obs_site = obs_site/Ns
    observables(1, head, bin) = observables(1,head, bin) + obs_site

      !! obs2: phonon kinetic energy
    obs_site = 0d0
    obs_site = obs_site  - 0.5d0 * M * sum((boson_field(:,time) - boson_field(:,time+1))**2)/(delt)**2
    obs_site = obs_site/Ns  + 1d0/(2d0*delt)
    observables(2, head, bin) = observables(2,head, bin) + obs_site

    do c1 = 1, Lat%N_cell
      pc1 => p_cells(c1)
      s1 = pc1%sites(1)%id
      !! obs3: electron kinetic energy along x direction
      obs_site = 0d0
      factor = (-1d0)**(1+(s1-1)/La)
      s1_next = s1 + 2
      if (mod(s1_next, La) == 1) s1_next = s1_next - La
      obs_site = obs_site + factor * (hop) * 2 * real(g_h(s1_next,s1) + g_h(s1,s1_next))
      observables(3,head, bin) = observables(3,head, bin) + obs_site/Ns

      !! obs4 : electron kinetic energy along y direction
      obs_site = 0d0
      factor = (-1d0)**(1+(s1-1)/La)
      s1_next = s1 + La
      if(s1_next > Ns) s1_next = s1_next - Ns
      s1_next = s1_next + 1
      if (mod(s1_next, La) == 1) s1_next = s1_next - La
      obs_site = obs_site + factor * (-hop)  * 2 * real(g_h(s1_next,s1) + g_h(s1,s1_next))
      observables(4,head, bin) = observables(4,head, bin) + obs_site/Ns

      !! obs5: electron potential energy
      obs_site = 0d0
      !obs_site = obs_site + 2 * real(ep_parameter * boson_field(s1,time) * (g_h(s1,s1)))
      !obs_site = obs_site -  real(ep_parameter**2 * g_h(s1,s1) * g_h(s1,s1)) - ep_parameter**2 * g_h(s1,s1)
      obs_site = obs_site  +  real(ep_parameter)**2 * g_h(s1,s1) * conjg(g_h(s1,s1))
      observables(5,head, bin) = observables(5,head, bin) + obs_site/Ns
      !! obs6: electron-phonon energy
      obs_site = 0d0
      obs_site = obs_site + 2 * real(ep_parameter * boson_field(s1,time) * (g_h(s1,s1)))
      observables(6,head, bin) = observables(6,head, bin) + obs_site/Ns
      !! obs7: electron total energy  
      obs_site = 0d0
      factor = (-1d0)**(1+(s1-1)/La)
      s1_next = s1 + 1
      if (mod(s1_next, La) == 1) s1_next = s1_next - La
      obs_site = obs_site + factor * (-hop) * 2 * real(g_h(s1_next,s1) + g_h(s1,s1_next))
      factor = (-1d0)**(1+(s1-1)/La)
      s1_next = s1 + La
      if(s1_next > Ns) s1_next = s1_next - Ns
      s1_next = s1_next + 1
      if (mod(s1_next, La) == 1) s1_next = s1_next - La
      obs_site = obs_site + factor * (-hop)  * 2 * real(g_h(s1_next,s1)+ g_h(s1,s1_next))
      if(M < 0.0001d0) then
        ! hubbard model
        obs_site = obs_site -  real(ep_parameter)**2 * g_h(s1,s1) * conjg(g_h(s1,s1))
      else
        !electron-phonon model
        obs_site = obs_site + 2 * real(ep_parameter * boson_field(s1,time) * (g_h(s1,s1)))
      end if
      observables(7,head, bin) = observables(7,head, bin) + obs_site/Ns

      !! obs8 : ph_X
      obs_site = 0d0
      obs_site = obs_site + (boson_field(s1,time))/Ns
      observables(8,head, bin) = observables(8,head, bin) + obs_site
      !! obs9: electron density
      obs_site = 0d0
      obs_site = obs_site + 2 * real(g_h(s1,s1))
      observables(9,head, bin) = observables(9,head, bin) + obs_site
      !! obs10: double occupancy
      obs_site = 0d0
      obs_site = obs_site + (g_h(s1,s1) * conjg(g_h(s1,s1)))
      observables(10,head, bin) = observables(10,head, bin) + obs_site/Ns
    end do !! for s1 and obs
    if(meas_corf) then
    do s1 = 1 , La
      do s2 = 1,corf_length
        corf_s2_list(s2,1) = modi(s1 + s2-1,La)
        corf_s2_list(s2,2) = corf_s2_list(s2,1) + La
      end do
      corf_s2_list(:,3:4) = corf_s2_list(:,1:2)
      corf_s1(1) = s1
      corf_s1(2) = s1
      corf_s1(3) = s1 + La
      corf_s1(4) = s1 + La
      do j = 1, 4 ! for A-A, A-B, B-A, B-B

            !! corf1: A-A 1-particle
        corf_offset = (j-1) * n_suit_corf
        corf_bin(:,1+corf_offset,bin) = corf_bin(:,1+corf_offset,bin) + g(corf_s2_list(:,j),corf_s1(j))/La

            !! corf2: A-A 2-particle
        corf_bin(:,2+corf_offset,bin) = corf_bin(:,2+corf_offset,bin) &
        & + g(corf_s2_list(:,j),corf_s1(j))*conjg(g(corf_s2_list(:,j),corf_s1(j)))/La

            !! corf3: nu_nu
        corf_bin(:,3+corf_offset,bin) = corf_bin(:,3+corf_offset,bin) + &
        & g_h(corf_s1(j),corf_s1(j))*cdiag(g_h(corf_s2_list(:,j),corf_s2_list(:,j))) /La
        corf_bin(:,3+corf_offset,bin) = corf_bin(:,3+corf_offset,bin) +&
        & g_h(corf_s1(j),corf_s2_list(:,j))* g(corf_s2_list(:,j),corf_s1(j)) /La

            !! corf4: nu_nd
        corf_bin(:,4+corf_offset,bin) = corf_bin(:,4+corf_offset,bin) + &
        & g_h(corf_s1(j),corf_s1(j))*conjg(cdiag(g_h(corf_s2_list(:,j),corf_s2_list(:,j))))/La

            !! corf5: den-den
        corf_bin(:,5+corf_offset,bin) = corf_bin(:,5+corf_offset,bin) +&
        & 4 * real(g_h(corf_s1(j),corf_s1(j))) * real(cdiag(g_h(corf_s2_list(:,j),corf_s2_list(:,j))))/La
        corf_bin(:,5+corf_offset,bin) = corf_bin(:,5+corf_offset,bin) + &
        & 2 * real(g_h(corf_s2_list(:,j),corf_s1(j))* g(corf_s1(j),corf_s2_list(:,j)))/La

            !! corf6: X-X
        corf_bin(:,6+corf_offset,bin) = corf_bin(:,6+corf_offset,bin) + boson_field(corf_s1(j),time)&
                                       &  * boson_field(corf_s2_list(:,j),time)/La

      end do ! for A-A, A-B, B-A, B-B
    end do !for s1 and corf
    end if ! FOR corf_measure

      !! END MEASUREMENT
    ! next loop
    if ((forward .and. time == ntime/2 + meas_number) .or. ((.not. forward) .and. time == ntime/2 - meas_number)) then
      ! average over all sites and all measurements in one loop
      observables(:,head, bin) = observables(:,head, bin)/(2*meas_number + 1)
      ! corf_bin between observables

      head = head + 1
      !print*,observables(head,bin,9)

      if (head > one_bin) then
        ! average the corf in one bin
        head = 1
        bin = bin + 1
        !if(record_middata) call mid_data_output(bin - 1)
      end if
    end if
    !2nd trotter
    if(second_trotter) then

    end if
  end subroutine

#ifdef MEASURE_TAU
  subroutine measure_along_tau()
    ! measure the observables that involve green funcitons with displaced time
    ! only enters here if forward
    ! if forward, time = ntime/2-ntau/2, G = G(time+1) = G(ntime/2-ntau/2+1)
    ! So G = G(ntime/2+1-ntau/2) which ends with exp(-delt*H(ntime/2+1-ntau/2)) at rightmost
    ! measured range is from ntime/2+1-ntau/2 to ntime/2+1+ntau/2
    ! e.g. ntime = 100, ntau = 20, G(41 : 61) with rightmost exp(-delt*H(41:61)),
    ! the invloved B(t + delt,t) is from (41,42) to (60,61), ntau = 20 in total
    ! we need the G(t,t),G_h(t,t) for t = 41,42,...,61,  G(t+1,t),G_h(t+1,t) for t = 41,42,...,60
    ! ngroup = 5, then ntau_block = 4, so the starting block is 9(t = 36:40), ending block is 13(t = 55:60)
    ! which means Q_string(9) is used for G(41,41), and also evolved from here
    implicit none
    integer :: tau, s1, s2, i,i_tau,ntau_block,bond_1,bond_2,s1_next,s2_next,j,b1,b2
    integer :: starting_block, ending_block, corf_site, bond_list(4,2)
    integer,allocatable :: block_list(:),tau_list(:)
    complex(8),allocatable ::  lmat(:,:), rmat(:,:)
    complex(kind=8) :: factor ! phase factor of the position (should be complex in general)
    complex(8) :: obs_site,  t_s1_c, t_s2_c
    real(8) ::  del_w, del_kx
    real(8) :: t_s1
    ! initialization
    ntau_block = ntau/ngroup
    starting_block = nblock/2 - ntau_block/2
    ending_block = nblock/2 + ntau_block/2
    !print*,'ntau_block:',ntau_block,'starting_block:',starting_block,'ending_block:',ending_block
    ! declaration
    if(.not. allocated(g_h_tau)) allocate (g_h_tau(Ns,Ns,ntau+1),g_tau(Ns,Ns,ntau+1))
    if(.not. allocated(g_h_dtau)) allocate (g_h_dtau(Ns,Ns,ntau),g_dtau(Ns,Ns,ntau))
    if(.not. allocated(g_h_tt)) allocate (g_h_tt(Ns,Ns,ntau+1),g_tt(Ns,Ns,ntau+1))
    if(.not. allocated(block_list)) allocate (block_list(ntau_block))
    if(.not. allocated(tau_list)) allocate (tau_list(ntau))
    do i = 1 , ntau
      ! tau_list = 41,42,...,60
      tau_list(i) = ntime/2 + i  - ntau/2
    end do

    ! initialize the PBP matrix and green function at tau = 1
    ! here g and g_h are at time = ntime/2+1-ntau/2, e.g. 41

    if(.not. allocated(lmat)) allocate (lmat(nelec,Ns),rmat(Ns,nelec))
    rmat = Q_string(:, :, starting_block) ! e.g.,starting block is 9
    lmat = transpose(conjg(Q_string(:, :, starting_block+1)))! e.g.,starting block+1 = 10

    g_h_tau(:,:,1) = g_h(:,:)
    g_tau(:,:,1) = g(:,:)
         !! calculate all g(tau,tau) and g(tau+1,tau)
    do i  = 1 , ntau
      ! i = 1:20 , tau = 41,42,...,60
      ! calculate G(tau + 1, tau) and G_h(tau + 1, tau)
      tau = tau_list(i)
      g_h_dtau(:,:,i) = g_h_tau(:,:,i)
      g_dtau(:,:,i) = g_tau(:,:,i)
      call evolve_tau(tau, g_dtau(:, :,i), ns, left = .true.,right = .false.)
      call evolve_tau(tau, g_h_dtau(:, :,i), ns, left = .false.,right = .true.)

      ! calculate G(tau + 1, tau+1) and G_h(tau + 1, tau+1)
      g_h_tau(:,:,i+1) = g_h_dtau(:,:,i)
      g_tau(:,:,i+1) = g_dtau(:,:,i)
      call evolve_tau(tau, g_tau(:, :,i+1), ns, left = .false.,right = .true.)
      call evolve_tau(tau, g_h_tau(:, :,i+1), ns, left = .true.,right = .false.)

            !! debug
      !call left_evolve(tau, g(:, :), ns, two_way = .true.)
      !call left_evolve_K(g(:, :), ns,two_way =.true.,half = .false.)
      !g_h(:, :) = idx(ns) - g(:, :)
      !print*,'fast_error:',maxval(abs(g_h(:,:) - g_h_tau(:,:,i+1)))
      ! stop
      ! numerical instability
      if(mod(i,ngroup) ==0) then
        ! tau = 45,50,55,60
        call get_g_scratch_tau(tau+1, rmat,g_h_tau(:,:,i+1))
        g_tau(:, :,i+1) = idx_c(Ns) - g_h_tau(:,:,i+1)
        !debug = .true.
        !call get_g_scratch_T0(tau+1, .true.,1)
        !debug = .false.
        !print*,'scratch_error:',maxval(abs(g_h(:,:) - g_h_tau(:,:,i+1)))

      end if
    end do

    ! initialize the g_tt and g_h_tt
    ! g_tt and g_h_tt are the displaced time green function at (1,tau), tau from 1 to ntau + 1
    ! when tau = 1, g_tt = g_tau(1,1), g_h_tt = g_h_tau(1,1)
    g_tt = 0d0
    g_h_tt = 0d0
    g_tt(:,:,1) = g_tau(:,:,1)
    g_h_tt(:,:,1) = g_h_tau(:,:,1)
    ! calculate g_h_tt and g_tt

    do i  = 1 , ntau
      g_tt(:,:,i+1) = matmul(g_dtau(:,:,i), g_tt(:,:,i))
      g_h_tt(:,:,i+1) = matmul(g_h_tt(:,:,i), g_h_dtau(:,:,i))
    end do

         !! let's measure

    do i_tau = 1 , ntau+1
      ! tau = 41, ... , 61
      !print*,real(g_tt(1,1,1,:))

      do s1 = 1, La
        do s2 = 1, La

          bond_list(:,1) = (/s1,s1,s1+La,s1+La/)
          bond_list(:,2) = (/s2,s2+La,s2,s2+La/)
          do j = 1, 4
            obs_site = 0d0
            bond_1 = bond_list(j,1)
            bond_2 = bond_list(j,2)
            obs_site = obs_site + cur_tau(i_tau,bond_1,bond_2,phase = .true.)
                     !! ctau 1: JJw_corf
            !factor = complex(0d0, del_w * (i_tau - 1))
            ctau_bin(i_tau,j,bin_tau) = ctau_bin(i_tau,j,bin_tau) +  delt * obs_site/La

                     !! corf 8: JJ_corf
            corf_site = s2 - s1
            if (corf_site < 0) corf_site = corf_site + La
            corf_site = corf_site + 1
            corf_bin(corf_site,j*8,bin_tau) = corf_bin(corf_site,j*8,bin_tau) + delt * obs_site/La
          end do
          !print*,'s2,s1:',s2,s1
          !print*,'obs_site:',obs_site

          ! JTx and JTw

          bond_list(:,1) = (/s1,s1+La,s1+2*La,s1+3*La/)
          bond_list(:,2) = (/s2,s2+La,s2+2*La,s2+3*La/)
          obs_site = 0d0
          do b1 = 1, 4
            do b2 = 1, 4
              bond_1 = bond_list(b1,1)
              bond_2 = bond_list(b2,2)
              obs_site = obs_site + cur_tau(i_tau,bond_1,bond_2,phase = .true.)
            end do
          end do
          ! JTx
          ctau_bin(i_tau,4 * n_suit_ctau + 1,bin_tau) = ctau_bin(i_tau,4 * n_suit_ctau + 1,bin_tau) +  delt * obs_site/La
          ! JTw
          corf_site = s2 - s1

          if (corf_site < 0) corf_site = corf_site + La
          corf_site = corf_site + 1
          corf_bin(corf_site,4 * n_suit_corf + 1,bin_tau) =&
          & corf_bin(corf_site,4 * n_suit_corf + 1,bin_tau) + delt * obs_site/La

        end do
      end do

      !print*,'ctau_bin:',ctau_bin(i_tau,5,bin_tau)
    end do
    !print*,'sum of ctau_bin:',sum(ctau_bin(:,5,bin_tau))
    !print*,'JJt',real(ctau_bin(:,5,bin_tau))
    !stop
  !! END MEASUREMENT_tau
    ! next loop
    ! corf_bin between observables

    head_tau = head_tau + 1
    !print*,observables(head,bin,9)

    if (head_tau > one_bin_tau) then
      head_tau = 1
      bin_tau = bin_tau + 1
      !if(record_middata) call mid_data_output(bin - 1)
    end if

  end subroutine measure_along_tau

  complex(8) FUNCTION cur_tau(time,bond_1,bond_2,phase) result(re)
    IMPLICIT NONE
    ! time is the time-index of the green function, from 1 to ntau+1
    ! s1,s2 are the index of the site, s1_next is the next site of s1
    ! s2_next is the next site of s2
    ! bond_1 and bond_2 are the index of the bond
    ! phase is logical, if .true. , the phase factor is included
    integer,intent(in) ::  bond_1, bond_2,time
    integer :: s1, s2, s1_next, s2_next
    logical :: phase
    integer :: flv
    complex(8) :: t_s1, t_s2
    re = 0d0
    s1 = bond_list(bond_1)%site_1
    s2 = bond_list(bond_2)%site_1
    s1_next = bond_list(bond_1)%site_2
    s2_next = bond_list(bond_2)%site_2

    t_s1 = get_bond_coupling(bond_1, time,phase = phase)
    t_s2 = get_bond_coupling(bond_2, time,phase = phase)
    do flv = 1 , nflv
      re = re + 4*real(( G_h_tau(s1_next,s1,time)-g_h_tau(s1,s1_next,time))) &
      & * real((G_h_tau(s2_next,s2,1)-g_h_tau(s2,s2_next,1))) + &
      & real(2*G_h_tt(s2_next,s1,time)*G_tt(s2,s1_next,time) + 2 * G_h_tt(s2,s1_next,time)*G_tt(s1,s2_next,time)-&
      & 2*G_h_tt(s2,s1,time)*G_tt(s1_next,s2_next,time) - 2*G_h_tt(s2_next,s1_next,time)*G_tt(s1,s2,time))
      re = -re * t_s1*t_s2

    end do

  end function cur_tau

  Subroutine evolve_tau(time, mat, n, left, right)
    ! evolve the matrix in the time+1 direction, special case for G(tau+1,tau)
    ! basically it is the same as left_evolve, but with a special case for G(tau+1,tau)
    implicit none
    integer,intent(in) :: n,time
    logical,intent(in) :: left,right
    complex(8),intent(inout) :: mat(ns,ns)
    integer :: site,pla,ph_i
    complex(8) :: expV

    if(left .and. right) then
      call left_evolve(time, mat, n, two_way=.true.)
      return
    end if
    if(left) then
      call left_evolve(time, mat, n, two_way=.false.)
    end if
    if(right) then
      ! special case for the G(tau+1,tau) calculation
      do site = 1, Ns
        call get_expV(expV, site, time)
        mat(:, site) = mat(:, site)/ expV
      end do
      if(.not.checkerboard) then
        mat = matmul(mat,expK_inv)
      else ! checherboard
        do ph_i = 1 , cb_field_number
        do pla = 1, n_cb_pla
          if(cb_pla_list(pla)%field_index .ne. ph_i) cycle
          mat(:,cb_pla_list(pla)%site_list) = matmul(mat(:,cb_pla_list(pla)%site_list),expK_inv)
        end do
        end do
      end if
    end if

  end subroutine

#endif

  subroutine data_analyse()
    ! analyse the data in the bins, actually the result of err part is not used

    IMPLICIT NONE
    integer :: i,corf
    if (.not. allocated(aver_bin)) allocate (aver_bin(nbin_per_core, n_obs))
    if (.not. allocated(obs_aver)) allocate (obs_aver(n_obs))
    if (.not. allocated(obs_err)) allocate (obs_err(n_obs))
    if (.not. allocated(corf_aver)) allocate (corf_aver(corf_length,n_corf), corf_err(corf_length,n_corf))
    if (.not. allocated(corfk_aver)) allocate (corfk_aver(corf_length,n_corf), corfk_err(corf_length,n_corf))
    aver_bin = 0d0
    obs_aver = 0d0
    obs_err = 0d0
    corf_aver = 0d0
    corf_err = 0d0
    corfk_aver = 0d0
    corfk_err = 0d0

    ! observables
    do obs = 1, n_obs
      aver_bin(:, obs) = sum(observables(obs,:, :), 1)/one_bin
      ! print*,'aver_bin :', aver_bin
      obs_aver(obs) = sum(aver_bin(:, obs))/nbin_per_core
      obs_err(obs) = sqrt(sum(aver_bin(:, obs)**2 - obs_aver(obs)**2)/(nbin_per_core - 1)/nbin_per_core)
    end do
         !! special obs
    aver_bin(:,11) = real(global_accept)/real(global_accept+global_reject)
    ! corf_bin
    if(meas_corf) then
      ! average the corf in one bin
      corf_bin = corf_bin/(2*meas_number + 1)/one_bin
      ! special normalization for the corf that takes one_bin_tau
      do i = 1 , 4
        corf_bin(:,i*8,:) = corf_bin(:,i*8,:)*(real(one_bin) * (2 *meas_number+1) /one_bin_tau)
      end do
      corf_bin(:,4 * n_suit_corf + 1,:) = &
      & corf_bin(:,4 * n_suit_corf + 1,:)*(real(one_bin) * (2 *meas_number+1) /one_bin_tau)
      ! calculate the corfk from corf
      do corf = 1, n_corf
        do bin = 1, nbin_per_core
          corfk_bin(:, corf, bin) = abs(matmul(corfk_phase,corf_bin(:, corf, bin)))
        end do
      end do
      ! transfer ctau into frequency space
      ctau_bin = ctau_bin / (one_bin_tau)
      ! special care for the crazy bins
      do bin = 1,nbin_per_core
        if(abs(ctau_bin(1,5,bin))>10) then
          ctau_bin(:,5,bin) = ctau_bin(:,5,modi(bin-1,nbin_per_core))
        end if
      end do
      !print*,'ctau_bin:',abs(ctau_bin(:,5,1))
      do corf = 1, n_ctau
        do bin = 1, nbin_per_core
          cw_bin(:, corf, bin) = real((matmul(ctau_phase,ctau_bin(:, corf, bin))))
        end do
      end do
      !print*,'ctau_bin:',abs(ctau_bin(:,5,1))

      ! average the corf and corfk and ctau in one core
      do corf = 1, n_corf
        do i = 1, corf_length
          corf_aver(i, corf) = sum(corf_bin(i, corf, :))/nbin_per_core
          corf_err(i, corf) =  sqrt(sum(corf_bin(i, corf, :)**2 - corf_aver(i, corf)**2)/(nbin_per_core - 1)/nbin_per_core)
          corfk_aver(i, corf) = sum(corfk_bin(i, corf, :))/nbin_per_core
          corfk_err(i, corf) = sqrt(sum(corfk_bin(i, corf, :)**2 - corfk_aver(i, corf)**2)/(nbin_per_core - 1)/nbin_per_core)
        end do
      end do

      do corf = 1, n_ctau
        do i = 1, ntau + 1
          ctau_aver(i, corf) = sum(ctau_bin(i, corf, :))/nbin_per_core
          ctau_err(i, corf) = sqrt(sum(ctau_bin(i, corf, :)**2 - ctau_aver(i, corf)**2)/(nbin_per_core - 1)/nbin_per_core)
          cw_aver(i, corf) = sum(cw_bin(i, corf, :))/nbin_per_core
          cw_err(i, corf) = sqrt(sum(cw_bin(i, corf, :)**2 - cw_aver(i, corf)**2)/(nbin_per_core - 1)/nbin_per_core)
        end do
      end do

    end if

    !print*,'obs_aver: ', obs_aver
    !print*,'obs_err: ', obs_err
  end subroutine data_analyse

  subroutine mid_data_output(tar_bin)
    IMPLICIT NONE
    integer :: flag = 0 , tar_bin
    character(len=20) :: ci1
    real(8) :: average_bin
    write (ci1, '(1i4)') modi(myid, MPI_one_block)

    do obs = 1, n_obs
      average_bin = sum(observables(obs,:,tar_bin))/one_bin
      flag = 0
      ci = obs_name(obs)
#ifdef MPI
      write (ci2, '(1i4)') MPI_block
      str = 'm_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'_core'//trim(adjustl(ci1))//'.dat'
#else
      str = 'm_'//trim(adjustl(ci))//'.dat'
#endif
      str = trim(adjustl(file_name))//'middata/'//trim(adjustl(str))
      inquire (file=str, exist=file_exist)
101   flag = flag + 1
      if(flag > 100) cycle
      if (.not. file_exist) then
        open (80 + obs, file=str, status='replace')
        close (80 + obs)
      end if
      !write(80+i,'(1f18.6)') 1d0/beta
      open (80 + obs, file=str, status='old', position='append',ERR = 101)
      write(80+obs,'(f18.6)') average_bin
      close (80 + obs)
    end do
  end subroutine

  subroutine aver_data_output()
    ! output the average data for each bin, core
    ! the obs data is going to be output in the form of o_obsname_b$block$.dat
    ! the corf data is going to be output in the form of c_corfname_b$block$c$core$.dat
    IMPLICIT NONE
    integer :: flag = 0,corf,i
    character(len=20) :: ci1
    ! output obs_aver and err
    do obs = 1, n_obs
      flag = 0
      ci = obs_name(obs)
#ifdef MPI
      write (ci2, '(1i4)') MPI_block
      str = 'o_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'.dat'
#else
      str = 'o_'//trim(adjustl(ci))//'.dat'
#endif
      str = trim(adjustl(file_name))//'data/'//trim(adjustl(str))
105   inquire (file=str, exist=file_exist)
      flag = flag + 1
      if(flag > 10) return
      if (.not. file_exist) then
        open (80 + obs, file=str, status='replace')
        close (80 + obs)
      end if
      open (80 + obs, file=str, status='old', position='append',ERR = 105)
      WRITE (ci, '(1i4)') nbin_per_core
      write (80 + obs, '('//trim(adjustl(ci))//'f18.12)',err = 106) real(aver_bin(:, obs))
106   close (80 + obs)
      if (myid < MPI_ONE_block .and. obs == 1) then
        !print *, myid, real(aver_bin(:, obs))
      end if
    end do

    ! output corf_aver and err
    if (meas_corf) then
    do corf = 1, n_corf
      flag = 0
      ci = corf_name(corf)
#ifdef MPI
      write (ci2, '(1i4)') MPI_block
      write (ci1, '(1i4)') mod(myid, MPI_one_block)
      str = 'c_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'c'//trim(adjustl(ci1))//'.dat'
#else
      str = 'c_'//trim(adjustl(ci))//'.dat'
#endif
      str = trim(adjustl(file_name))//'data/'//trim(adjustl(str))
107   inquire (file=str, exist=file_exist)
      flag = flag + 1
      if(flag > 10) return
      if (.not. file_exist) then
        open (80 + obs, file=str, status='replace')
        close (80 + obs)
      end if
      open (80 + obs, file=str, status='old', position='append',ERR = 107)
      WRITE (ci, '(1i4)') nbin_per_core
      write (ci1, '(1i4)') corf_length
      do i = 1, corf_length
        write (80 + obs, '('//trim(adjustl(ci))//'f18.12)',err = 108) real(corf_bin(i, corf, :))
      end do
108   close (80 + obs)
    end do

         !! output corfk_aver and err
    do corf = 1, n_corf
      flag = 0
      ci = corf_name(corf)
#ifdef MPI
      write (ci2, '(1i4)') MPI_block
      write (ci1, '(1i4)') mod(myid, MPI_one_block)
      str = 'ck_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'c'//trim(adjustl(ci1))//'.dat'
#else
      str = 'ck_'//trim(adjustl(ci))//'.dat'
#endif
      str = trim(adjustl(file_name))//'data/'//trim(adjustl(str))
109   inquire (file=str, exist=file_exist)
      flag = flag + 1
      if(flag > 10) return
      if (.not. file_exist) then
        open (80 + obs, file=str, status='replace')
        close (80 + obs)
      end if
      open (80 + obs, file=str, status='old', position='append',ERR = 109)
      WRITE (ci, '(1i4)') nbin_per_core
      write (ci1, '(1i4)') corf_length
      do i = 1, corf_length
        write (80 + obs, '('//trim(adjustl(ci))//'f18.12)',err = 110) abs(corfk_bin(i, corf, :))
      end do
110   close (80 + obs)
    end do
    end if

         !! output ctau_aver and err
    if ( meas_corf) then
    do corf = 1, n_ctau
      flag = 0
      ci = ctau_name(corf)
#ifdef MPI
      write (ci2, '(1i4)') MPI_block
      write (ci1, '(1i4)') mod(myid, MPI_one_block)
      str = 'ct_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'c'//trim(adjustl(ci1))//'.dat'
#else
      str = 'ct_'//trim(adjustl(ci))//'.dat'
#endif
      str = trim(adjustl(file_name))//'data/'//trim(adjustl(str))
      flag = flag + 1
      if(flag > 10) return
111   inquire (file=str, exist=file_exist)
      if (.not. file_exist) then
        open (80 + obs, file=str, status='replace')
        close (80 + obs)
      end if
      open (80 + obs, file=str, status='old', position='append',ERR = 111)
      WRITE (ci, '(1i4)') nbin_per_core
      write (ci1, '(1i4)') ntau + 1
      do i = 1, ntau + 1
        write (80 + obs, '('//trim(adjustl(ci))//'f18.12)',err = 112) real(ctau_bin(i, corf, :))
      end do
112   close (80 + obs)
    end do
    end if

         !! output cw_aver and err
    if ( meas_corf) then
    do corf = 1, n_ctau
      flag = 0
      ci = ctau_name(corf)
#ifdef MPI
      write (ci2, '(1i4)') MPI_block
      write (ci1, '(1i4)') mod(myid, MPI_one_block)
      str = 'cw_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'c'//trim(adjustl(ci1))//'.dat'
#else
      str = 'cw_'//trim(adjustl(ci))//'.dat'
#endif
      str = trim(adjustl(file_name))//'data/'//trim(adjustl(str))
      flag = flag + 1
      if(flag > 10) return
113   inquire (file=str, exist=file_exist)
      if (.not. file_exist) then
        open (80 + obs, file=str, status='replace')
        close (80 + obs)
      end if
      open (80 + obs, file=str, status='old', position='append',ERR = 113)
      WRITE (ci, '(1i4)') nbin_per_core
      write (ci1, '(1i4)') ntau + 1
      do i = 1, ntau + 1
        write (80 + obs, '('//trim(adjustl(ci))//'f18.12)',err = 114) real(cw_bin(i, corf, :))
      end do
114   close (80 + obs)
    end do
    end if

  end subroutine aver_data_output

  function cdiag(A) result(re)
    complex(8), intent(in) :: A(:,:)
    complex(8):: re(size(A,1))
    integer :: i
    do i = 1, size(A,1)
      re(i) = A(i,i)
    end do
  end function

end module
