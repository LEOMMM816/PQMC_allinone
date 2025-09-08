program Pqmc_main

  use dqmc_measurements
  use update
  !use mpi

  IMPLICIT NONE
#ifdef MPI
  include 'mpif.h'
#endif
! parameters declare

  integer :: loop, time, time_count, iflv, i_phonon_field, i_site,i_bond,i_pf
  integer :: i_ml, i_ml2
  integer(4) :: count1, count2, count_rate
  integer :: iteration  ! # of Monte Carlo loops
  logical :: forward = .true.
  character*10 :: ci_myid

!> mpi
#ifdef MPI
  call MPI_INIT(ierr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
  CALL MPI_COMM_rank(MPI_COMM_WORLD, myid, ierr)
  MPI_block = myid/MPI_one_block + 1
  if(numprocs .NE. mpi_nblock*MPI_one_block) then
    print*,'numprocs:',numprocs,'mpi_nblock:',mpi_nblock,'MPI_one_block:',MPI_one_block
    stop 'numprocs!=mpi_nblock*MPI_one_block'
  end if
  write (ci_myid, '(i4.4)') myid
  nml_file = 'namelists/input'//trim(adjustl(ci_myid))//'.nml'
  if ( fixedseeds ) then
    call init_rng_by_seed(122)
    print *, 'seed:', ranseeds
  else
    call init_rng(myid + 1)
    print *, 'myid:',myid,'seed:', ranseeds
  end if
#endif

! run time
  call system_clock(count1, count_rate)

  ! lattice initialize
  print *, 'setting lattice...'
  call set_lattice_main()
  print *, 'setting lattice ends.'
  call init()
  print*,"setting boson field..."
  call init_boson_field()
  print*,"setting bf ends."
  ! phonon-field initialize
  print *, 'set phonon field...'
  call set_pf_main() 
  call readin_slater_pf()
  print *, 'set phonon field ends.'
#ifdef MPI
  call init_MPI(MPI_block)
#endif
  ! slater initialization
  call set_k_slater()
  call set_nelec_occupied()
  ! initialize for measurement
  call Meas_sys%init_meas()
#ifdef MPI 
  print *, 'myid:',myid, '/', &
  &'lattice_diemension:',Lat%dim, '/',&
  &'filling:',filling,'nelec:', nelec,'/', &
  &'beta,ntime,lambda,nblock', beta, ntime, ep_parameter, nblock
#else
  print *, 'nelec:', nelec
  print *, 'beta,ntime,nblock', beta, ntime, nblock
#endif
!simulation begins here
  call init_MC()
  
  do loop = 1, iteration ! Monte Carlo steps ()

    call middle_check()

    forward = .not. forward ! the first loop is going backwards

    ! global update
    if (.not. forward .and. mod(loop + 1, global_update_loop) == 0 .and. global_update) then
      updated = .false.

      global_update_flip = .not. global_update_flip
      global_update_shift = .not. global_update_shift
      global_update_exchange = .false.
      call update_global()

      global_update_exchange = .true.
      call update_global()

      if(updated) call init_g_T0()

    end if
    ! update every phonon field
    do time_count = 1, ntime

      !fast update green function along time direction (backward) ! at time to get g(time)
      if (.not. forward) then
        time = ntime - time_count + 1
        ! G(time+1) -> B^(-1) G(time+1) B = G(time),  where B = exp(-delt*H(time))
        do i_pf =  n_phonon_field,1,-1
          call right_evolve(time, g, ns,i_pf,two_way = .true.)
          g_h = idx(ns) - g
          if (local_update .and. ppf_list(i_pf)%V_exist) then
            ! update G(time) where the exp(-delt*V(time)) is to the right most
            call update_phonon_field(ppf_list(i_pf), time)
          end if

        end do ! i_pf
      end if ! backward

      !fast update green function along time direction (forward) ! at time to get g(time+1)
      if (forward) then
        time = time_count
        ! G(time) -> B G(time) B^(-1) = G(time+1),  where B = exp(-delt*H(time))
        do i_pf = 1, n_phonon_field
          if (local_update .and. ppf_list(i_pf)%V_exist) then
            ! update G(time) where the exp(-delt*V(time)) is to the right most
            call update_phonon_field(ppf_list(i_pf), time)
          end if
          call left_evolve(time, g, ns,i_pf,two_way = .true.)
          g_h = idx(ns) - g
        end do
      end if ! forward

      ! calculate green function from scratch if ngroup time slices have passed
      if (mod(time_count, ngroup) == 0) then

        if (forward) then
          !print *, 'forward time:', time
          call get_g_scratch_T0(time + 1, forward, loop) ! at time to get g(time+1)
        else
          call get_g_scratch_T0(time, forward, loop) ! at time to get g(time)
        end if

      end if ! end update g from scratch

! Measurement begins

      if (loop > warmup .and. mod(loop - warmup, meas_interval) == 0) then
        if(time >= ntime/2 - meas_number .and. time <= ntime/2 + meas_number) then
          ! take measurement every meas_interval MCS
          call Meas_sys%begin_measure( forward,time)
        end if
      end if

      if (loop > warmup .and. mod(loop - warmup, meas_interval_tau) == 0) then
        if(time_count == ntime/2-ntau/2 .and. forward) then
          ! if forward, time = ntime/2-ntau/2, G = G(time+1) = G(ntime/2+1-ntau/2)
          !  call measure_along_tau()
        end if
      end if ! for measurement

    end do ! for time in 1 :ntime

  end do ! for MCS

! data analysis
  call Meas_sys%data_analyse()
! wrint bin-averaged-per-core data into file
  !call Meas_sys%aver_data_output()
! write (mpi-block/core) averaged data into file

   !if(myid == 0) call MPI_output_final()
!write detailed data into file
 ! if (record_ph_field) call output_phonon_field()
! output info_sheet
#ifdef MPI
  if (mod(myid,MPI_one_block) == 0) call mpi_output_inforsheet()
#else
  call mpi_output_inforsheet()
#endif
! output error information
  print *, 'largest error in', myid, 'process is', err_fast_max
! run time
  call system_clock(count2)
  print *, 'RUN TIME(s): ', real(count2 - count1)/real(count_rate)

! mpi ends
#ifdef MPI
  call MPI_Finalize(ierr)
#endif

contains
  subroutine init_MC()
    implicit none
    if (.not. allocated(g)) then
    do i_pf = 1 , n_phonon_field
      call init_expKV(ppf_list(i_pf))
    end do
    call init_g_T0()
    print*,'init g_T0 finishes'
    end if

  iteration = (nbin_per_core*nmeas_per_bin)*meas_interval + warmup
  forward = .true.
    
  end subroutine

  subroutine middle_check()
    implicit none
    if (mod(loop, print_loop) == 0 .and. mod(myid, MPI_one_block) == 0) then
      inquire (file='check_middle.dat', exist=file_exist)

      if (.not. file_exist) then
        open (110, file='check_middle.dat', status='replace')
        close (110)
      end if
      open (110, file='check_middle.dat', status='old', position='append')
      WRITE (110, *) 'loop',loop
      WRITE (110, *) 'myid',myid
      write (110, *) 'Size:', Lat%dlength
      if(ML_update) then
        write (110, *) 'ML_accept :', real(ML_accept)/(ML_accept + ML_reject)
        write (110, *) 'ML_accept_ratio :', ML_accept_ratio / ML_accept
        write (110, *) 'ML_t_accept :', real(ML_t_accept)/(ML_t_accept + ML_t_reject)
        write (110, *) 'ML_t_accept_ratio :', ML_t_accept_ratio / ML_t_accept
        write (110, *) 'ML_s_accept :', real(ML_s_accept)/(ML_s_accept + ML_s_reject)
        write (110, *) 'ML_s_accept_ratio :', ML_s_accept_ratio / ML_s_accept
        ! change update_distance if needed
        if ( ml_s_accept_ratio/(ML_s_accept + ML_s_reject) < 0.1) then
          ml_distance_ratio = ml_distance_ratio * 0.5d0
        end if
        ! write (110, *) 'ML_total :', (ML_accept + ML_reject)
        write (110, *) 'wolff_accept :', wolff_accept
        write (110, *) 'wolff_total :', (wolff_accept + wolff_reject)
        !write (110, *) 'wolff_time_accept :', wolff_time_accept
      end if
      if(global_update) then
        write (110, *) 'global_update accept :', global_accept
        write (110, *) 'global_update total :', (global_accept + global_reject)
      end if
      !WRITE(110,'(1f18.6)')
      write (110, *) 'local_update accept ratio:', real(positive_accept)/(positive_accept + negative_accept)
      !WRITE(110,'(1f18.6)') real(positive_accept)/(positive_accept+negative_accept)
      write (110, *) 'error_fast:', err_fast_max
      !WRITE(110,'(1f18.6)') err_fast_max
      !write(110,*)'delta_energy_E,delta_energy_K,delta_energy_P:',delta_energy_E,delta_energy_K,delta_energy_P
      !WRITE(110,'(3f18.6)') delta_energy_E,delta_energy_K,delta_energy_P
      write(110, *) 'max_pair:',max_pair
      call system_clock(count2)
      WRITE (110, *) 'time spent', real(count2 - count1)/real(count_rate)
      write (110, *) ''
      close (110)
      positive_accept = 0
      negative_accept = 0
      delta_energy_E = 0d0
      delta_energy_K = 0d0
      delta_energy_P = 0d0
      ml_s_accept = 0
      ml_s_reject = 0
      ml_t_accept = 0
      ml_t_reject = 0
      ml_accept = 0
      ml_reject = 0
      ml_s_accept_ratio = 0d0
      ml_t_accept_ratio = 0d0
      ml_accept_ratio = 0d0

      !wolff_accept = 0
      !wolff_reject = 0

    end if
  end subroutine
! ////////////////////////////////////////////////////////////////////////////////////////

  SUBROUTINE set_k_slater()
    implicit none
    integer :: i_pla,i_pf
    real(8) ::  factor
    
    if (.not. allocated(K_slater)) ALLOCATE (K_slater(Ns, Ns))
    K_slater = 0
    ! build K_slater matrix from typed slater_pf
    do i_pla = 1, slater_pf%n_plaquette
      K_slater(slater_pf%p_data%pla_site_list(:,i_pla), slater_pf%p_data%pla_site_list(:,i_pla))&
      &  =  slater_pf%K_coe * slater_pf%Kmatrix
    end do

    
      allocate(K_mat(Ns,Ns),expK(Ns,Ns),expK_half(Ns,Ns),expK_inv(Ns,Ns),expK_inv_half(Ns,Ns))
      K_mat = 0d0
      do i_pf = 1, n_phonon_field
        if (ppf_list(i_pf)%K_exist) then
          do i_pla = 1, ppf_list(i_pf)%n_plaquette
          K_mat(ppf_list(i_pf)%p_data%pla_site_list(:,i_pla), ppf_list(i_pf)%p_data%pla_site_list(:,i_pla))&
      &  =  K_mat(ppf_list(i_pf)%p_data%pla_site_list(:,i_pla), ppf_list(i_pf)%p_data%pla_site_list(:,i_pla)) + &
      & ppf_list(i_pf)%K_coe * ppf_list(i_pf)%Kmatrix
      end do
        end if
      end do
      
      expK = -delt * K_mat
      call expm(expk, Ns)
      expK_half = -delt/2d0 * K_mat
      call expm(expK_half, Ns)
      expK_inv = delt * K_mat
      call expm(expK_inv, Ns)
      expK_inv_half = delt/2d0 * K_mat
      call expm(expK_inv_half, Ns)
    
    
  END SUBROUTINE

  subroutine set_nelec_occupied()

    IMPLICIT NONE
    INTEGER :: flv
    REAL(8) eval(Ns), eval_comple(Ns)

    !IF(abs(dot_product(expk_half(1,:,1),inv_expk_half(:,1,1))-1d0)>1d-6)THEN
    !   PRINT*,'expk_half has not been correctly set. It is required for 2nd-order Trotter.'
    !   CALL exit(0)
    !END IF
    IF (.not. allocated(slater)) ALLOCATE (slater(Ns, nelec))
    IF (.not. allocated(slater_Q)) ALLOCATE (slater_Q(Ns, nelec))
    IF (.not. allocated(slater_D)) ALLOCATE (slater_D(nelec))
    IF (.not. allocated(R_string)) ALLOCATE (R_string(nelec, nelec))
    !open(39,file='eval.dat')

    CALL eigen(Ns, K_slater(:, :), eval)

    slater(:, 1:nelec) = K_slater(:, 1:nelec)
    slater_Q(:, 1:nelec) = slater(:, 1:nelec)
    CALL qdr(Ns, nelec, slater_Q(:, 1:nelec), R_string(1:nelec, 1:nelec), slater_D(1:nelec))

    !close(39)
  end subroutine

  subroutine init_boson_field()

    integer :: i_cell, l, time, k_index,ios,i_bf
    logical :: ising,continuous
    real(8) :: random_range, offset
    complex(8),allocatable :: phase_intra_cell(:)
    real(8),allocatable :: phase_k_inter_cell(:)
    complex(8) :: phase_factor
    character*4 :: ci3
    character*30 :: str = 'ph_field_.dat'
    namelist /bf_basics/ bf_sets,ising,continuous,import_ph_field
    namelist /bf_initial/ random_range, offset,phase_intra_cell,phase_k_inter_cell
    open(10, file = nml_file, status='old', position='rewind')
    read(10, nml=bf_basics)
    n_boson_field = bf_sets * Lat%N_cell
    if (.not. allocated(boson_field)) allocate (boson_field(bf_sets*Lat%N_cell,ntime))
    boson_field = 0d0
    do i_cell = 1, Lat%N_cell
      p_cells(i_cell)%bf_list = [(i_bf + (i_cell - 1) * bf_sets, i_bf = 1, bf_sets)]
    end do
    if (.not. import_ph_field) then ! set up ph field from new
      allocate(phase_intra_cell(bf_sets), phase_k_inter_cell(Lat%dim))
      read(10, nml=bf_initial, iostat=ios)
      if (ios /= 0) then
        print *, "Error reading namelist bf_initial"
        print*,'ios:',ios
        stop
      end if
      random_range = random_range * char_length
      offset = offset * char_length
      do l = 1, ntime
        !boson_field(:,:,l,flv,i_pf) = K
        do i_cell = 1, Lat%N_cell
          phase_factor = exp(complex(0d0,sum(phase_k_inter_cell * p_cells(i_cell)%rpos * PI)))
          do i_bf = 1, bf_sets
            boson_field(p_cells(i_cell)%bf_list, l) = real((random_range*ran_sysm() + (phase_factor*offset))*phase_intra_cell(i_bf))

          end do
        end do
      end do

      !boson_field = 0d0
      !k_max = 40
      boson_field(:,1) = end_field
      boson_field(:,ntime) = end_field

    else ! import phonon field

#ifdef MPI
      write (ci3, '(1i4)') 1
      str = 'ph_field_core'//trim(adjustl(ci3))//'.dat'
#endif
      write (ci3, '(1i4)') bf_sets * Lat%N_cell
      open (99, file=str, status='old', position='rewind')

      do time = 1, ntime
        read (99, '('//trim(adjustl(ci3))//'f18.6)') boson_field(:, time)
      end do
      close (99)
    end if

  end subroutine

end program Pqmc_main
