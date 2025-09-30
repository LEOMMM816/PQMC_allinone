
MODULE input

  IMPLICIT NONE

! debug varibles
  integer :: positive_accept = 0, negative_accept = 0, global_accept = 0, global_reject = 0
  real(8) :: delta_energy_P = 0, delta_energy_K = 0, delta_energy_E = 0
  integer :: i_debug
! Control varibles
  logical :: proj = .true. ! pqmc or ftqmc
  logical :: record_ph_field = .false.
  logical :: Metropolis = .true.
  logical :: check_acorr = .false.
  logical :: import_ph_field = .false.
  logical :: greatest_decent = .false.
  logical :: local_update = .true.
  logical :: Kspace_GU = .false.
  logical :: x_y_inv = .false.
  logical :: ML_weight = .false.
  logical :: ML_update = .false.
  logical :: TBC = .false.
  logical :: global_update = .false.
  logical :: global_update_shift = .false.
  logical :: global_update_flip = .true.
  logical :: global_update_exchange = .false.
  logical :: record_middata = .false.
  logical :: debug = .false.
  logical :: meas_corf = .true.
  logical :: fixedseeds = .false.
  logical :: second_trotter = .false.

! constant
  real(8),parameter :: e = 2.7182818284590452353602874713526624d0

! mpi parameter
  integer :: ierr, myid = 0, numprocs = 1, MPI_block = 1
  integer :: MPI_one_block = 1 !> number of cores in one mpi_block, = 1 if not mpi
  integer :: MPI_nblock = 1!> = numprocs/MPI_one_block
! MC parameter
  integer :: warmup = 2000, meas_interval = 1, meas_interval_tau = 20, meas_number = 0
  integer :: nbin_per_core = 20, nmeas_per_bin = 1000! # of bins and size of one bin
! measurement parameter
  integer :: n_suit_corf = 8 , n_suit_ctau = 1 ! number of suited correlation functions
  logical :: file_exist
  character*30, allocatable :: obs_name(:),corf_name(:),ctau_name(:)
  integer,allocatable :: corf_length(:)
  integer :: ntau = 80, one_bin_tau

! General ARGUMENTS
  character(len=100) :: nml_file = 'input.nml', output_addr = 'data/',output_format = 'f18.8',reading_guide = 'reading_guide.nml'
  character(len=200) :: task_name
  integer :: nflv = 1, ncopy = 1 ! number of flavours
  integer :: ntime, Ns ! ntime is number of time slices and Ns is # of sites
  integer :: nblock, ngroup = 5 !> nblock = ntime/ngroup + 2
  real(8) :: delt = 0.1, beta = 80d0,hop = 1d0
  integer :: print_loop = 20
  integer :: lattice_dimension
  real(8) :: prob !> accept probobility
  real(8) :: err_fast = 0d0, err_fast_max = 0.000001d0
  integer, allocatable :: field_type(:) ! field number
  real(8):: TBC_phase_number(2) = (/0d0,0d0/)
  real(8) :: filling = 1d0
  real(8) :: biased_phonon

! HS parameter

  real(kind=8) :: ep_parameter = 1d0, U = 0d0 ! Hubbard U and steps delt * U < 0.5

! proj
  real(8) :: disorder = 0d0
  real(kind=8) :: eta ! et a = exp(-delt*(miu - U/2))
  real(8) :: filling_factor
  integer :: nelec
  real(kind=8) R_nflv ! relative boltzmann weight for different particles
  complex(8), ALLOCATABLE :: slater(:, :), slater_Q(:, :), slater_D(:)
  complex(kind=8), allocatable :: K_slater(:, :)! kinetic energy ns*ns

! phonon parameter
  integer :: n_phonon_field = 4 !> number of phonon fields
  integer :: bf_sets = 2 ! number of boson fields
  integer :: n_boson_field ! number of decomposed phonon field
  real(8) :: D = 1d0,M = 1d0 !>stiffness constant and Mass
  real(8) :: char_length = 1.0 !> characteristic length of phonon field
  real(8) :: omega = 1d0!> phonon frequency
  real(8) :: end_field = 0d0
  real(8) :: jump_distance
  real(8) :: max_displacement !> = hop / ep_parameter
  integer :: n_local_update = 1 !> times of proposal  of new ph field at each location
  real(8), allocatable,target :: boson_field(:,:) !(bf_sets*N_cell,time)

!global update
  real(8) :: ln_cw = 0d0
  real(8) ::  global_update_distance = 1.0d0 !>
  integer :: n_global_update = 1!>number of bonds that change in one GUD
  integer :: global_update_loop = 2
  logical :: updated = .false.
!ML update
  integer :: ML_accept = 0, ML_reject = 0, ML_t_accept = 0, ML_t_reject = 0, ml_s_accept = 0, ml_s_reject = 0
  integer :: wolff_accept = 0, wolff_reject = 0 , wolff_time_accept = 0, wolff_time_reject = 0
  real(8) :: ML_accept_ratio = 0d0, ML_t_accept_ratio = 0d0, ml_s_accept_ratio = 0d0
  real(8) :: ml_distance_ratio = 0.75d0
  integer :: n_update_ml = 5, n_local_update_ml = 20
  integer :: n_k_ML = 40
! green function related
  complex(kind=8), allocatable :: K_mat(:, :),expK(:,:),expK_half(:,:),expK_inv(:,:), expK_inv_half(:,:)
  complex(kind=8), allocatable :: g(:, :), g_h(:, :) ! green function(ns,ns) & inv
  !complex(kind = 8), allocatable :: g_debug(:,:,:)
  complex(8), ALLOCATABLE ::  Q_string(:, :, :), D_string(:, :) ! (Ns),nelec,nblock
  complex(8), ALLOCATABLE :: R_string(:, :) ! auxilliary matrix to store R_matrix in qdr decomposition
  ! new things
  type :: pf_data_type
    integer,allocatable :: pla_site_list(:,:)! (pf%dim,n_plaquette) the site index in each plaquette for each phonon field
    integer,allocatable :: bf_list(:) !> n_plaquette, index of the coupled boson field
    complex(8),allocatable :: expKV(:,:,:,:) !> exp(-delt* K or V) for each pf, (ppf%dim,ppf%dim,n_plaquette,ntime)
    complex(8),allocatable :: expKV_inv(:,:,:,:) !> exp(delt* K) for each pf, (ppf%dim,ppf%dim,n_plaquette,ntime)
  end type pf_data_type
  type :: lat_type
    logical :: periodic
    integer :: dim,n_subsite_uc,N_cell,Ns
    integer, allocatable :: dlength(:)
    real(8), allocatable :: rlength(:) ! length of the lattice in each dimension
    real(8),allocatable :: tsl_rvec(:,:),subsite_rvec(:,:)
    ! tsl_rvec is the lattice vectors in spatial vectors
    ! subsite_rvec is the subsites' vectors in spatial vectors

  end type lat_type

  type :: pf_type
    integer :: id
    logical :: K_exist, V_exist
    integer :: dim
    integer :: n_plaquette
    complex(8) :: K_coe, V_coe !> coefficient of K and V in the pf Hamiltonian
    complex(8),allocatable :: Vmatrix(:,:)
    complex(8),allocatable :: Kmatrix(:,:)
    integer, allocatable :: pla_tsl_dvec_uc(:,:)
    !> translational vectors linking all unit cells hosting the first site,(lat%dim,lat%dim)
    integer, allocatable :: pla_offset_dvec_uc(:)
    !> offset translational vector (lat%dim)
    integer, allocatable :: pla_int_subsites(:,:)
    !>[d_vecs,sub_index] of sites in each pla, including the first site,(lat%dim+1,ppf%dim)
    type(pf_data_type),pointer :: p_data !> phonon field data
    logical,pointer :: debug
  end type pf_type

  type:: cell_type
    integer :: id! unique cell id
    integer, allocatable :: dpos(:) ! position of the unit cell in lattice vectors
    real(8), allocatable :: rpos(:) ! position in spatial vectors
    type(site_type), pointer :: sites(:) ! sites in the unit cell
    integer, allocatable :: bf_list(:) ! boson fields in this cell
  end type cell_type

  type :: site_type
    integer :: id,uc_id,subsite_id ! unique site id, unit cell id, subsite id
    integer,allocatable :: uc_dpos(:) ! position of the unit cell in lattice vectors
    real(8), allocatable :: uc_rpos(:),site_rpos(:) ! position in spatial vectors
    type(cell_type), pointer :: p_uc ! pointer to the cell this site belongs to
  end type site_type

  type(lat_type) :: Lat
  type(cell_type), allocatable,target :: typed_cells(:) ! cells in the lattice
  type(cell_type), pointer :: p_cells(:) ! pointer to cells in the lattice
  type(site_type), allocatable,target :: typed_sites(:)
  type(site_type), pointer :: p_sites(:)
  type(pf_type),allocatable,target :: pf_list(:)
  type(pf_type),pointer :: ppf_list(:)
  type(pf_type),target :: slater_pf
  type(pf_data_type),allocatable, target :: pf_data_list(:)
  type(pf_data_type),target :: slater_pf_data

contains

  SUBROUTINE init()
    implicit none
    ! parameters definition
    one_bin_tau = nmeas_per_bin/(meas_interval_tau/meas_interval)

    ntime = nint(beta/delt)
    Ns = Lat%Ns
    nblock = ntime/ngroup + 2
    if(M < 0.001d0) then
      omega = 100000000d0
      char_length = sqrt(1d0/(D*ep_parameter))
    else
      omega = sqrt(D/M)
      char_length = (1d0/(M*omega))
    end if
    jump_distance = (5d0) * (1d0/delt) * char_length
    filling_factor = filling
    nelec = nint(Ns*filling_factor)
    U = ep_parameter**2/D
  end subroutine

  subroutine init_MPI(mpi_block)
    integer, intent(in) :: mpi_block
    integer :: i_pf
    character(len=20) :: temp_string
    D = 1d0
    M = 0d0
    ep_parameter = sqrt(mpi_block + 2d0) ! hubbard U = ep^2
    filling = 10d0/40d0
    biased_phonon = 0d0
    call init()

    do i_pf = 1, n_phonon_field
      pf_list(i_pf)%K_coe = (-hop)
      pf_list(i_pf)%V_coe = ep_parameter
    end do
    slater_pf%K_coe = (-hop)
    return
  end subroutine init_MPI

!------------------------------------------------------!
!-------------------Output_info-------------------------!
!------------------------------------------------------!

  subroutine mpi_output_inforsheet()
    implicit none
      !! output all information of this simulation
      !! take the format of the La12Lb12Lt96dau... .out file in this folder
    character(len = 100) :: ci
    integer :: i_dim
    integer :: unit, iterations
    task_name = trim(output_addr)//'out_files/'
#ifdef MPI
    write(ci,'(a,1i4)') 'b',MPI_block
    task_name = trim(task_name)//trim(ci)
#endif
    task_name = trim(task_name)//'Len'
    do i_dim = 1, Lat%dim
      write(ci, '(I3)') Lat%dlength(i_dim)
      task_name = trim(task_name)//trim(ci)
      if(i_dim < Lat%dim) task_name = trim(task_name)//'*'
    end do
    write(ci,'(a,1i4,a,1f5.2,a,1i4)') 'Lt', ntime, 'dtau',delt,'Nup',nelec
    task_name = trim(task_name)//trim(ci)
    if(M > 0.0001d0) then
      ! phonon model
      write(ci, '(a,1f5.2)') 'lam', real(ep_parameter)
    else
      ! hubbard model
      write(ci, '(a,1f5.2)')'U', U
    end if
    task_name = trim(task_name)//trim(ci)//'.out'
    call remove_spaces(task_name,len(task_name))
    iterations = (nbin_per_core*nmeas_per_bin)*meas_interval + warmup
    open(newunit=unit, file=task_name, status='replace', action='write')

    ! 写入文件头
    write(unit, '(a)') repeat('*', 48)
    write(unit, '(a)') repeat(' ', 1) // repeat('*', 47)
    write(unit, '(a)') repeat(' ', 11) // trim(adjustl(nml_file)) // repeat(' ', 11) // '*'
    write(unit, '(a)') repeat('*', 48)
    write(unit, '(a)') repeat('*', 48)
    write(unit, '(a)') repeat(' ', 1) // 'Hamiltonian: H = H_k + \lambda * (n_up+n_down)'
    write(unit, '(a)') repeat(' ', 1) // 'Biased phonon used, H_ep = \lambda * X * (n_up + n_down - <n>)'
    write(unit, '(a)') ''

    ! put the monte carlo paramters
    write(unit, '(a)') repeat(' ', 1) // repeat('#', 36)
    write(unit, '(a)') repeat(' ', 1) //repeat('#', 5) // '  Monte Carlo parameter  ' // repeat('#', 5)
    write(unit, '(a)') repeat(' ', 1) // repeat('#', 36)
    write(unit, '(a, i15)') ' nbin  = ', nbin_per_core
    write(unit, '(a, i15)') ' measurements in one_bin', nmeas_per_bin
    write(unit, '(a, i15)') ' warmup  = ', warmup
    write(unit, '(a, i15)') ' meas_interval  = ', meas_interval
    write(unit, '(a, i15)') ' meas_timerange  = ', meas_number
    write(unit, '(a, i15)') ' total_iterations  = ', iterations
    write(unit, '(a, i15)') ' number of mpi_blocks  = ', MPI_nblock
    write(unit, '(a, i15)') ' number of cores per block  = ', MPI_one_block
    if(global_update) then
      write(unit, '(a)') ' global update is used'
    else
      write(unit, '(a)') ' global update is not used'
    end if
    write(unit, '(a)') ''
    ! 写入系统参数
    write(unit, '(a)') repeat('#', 36)
    write(unit, '(a)') repeat(' ', 1) // repeat('#', 7) // '  Physical parameter  ' // repeat('#', 7)
    write(unit, '(a)') repeat('#', 36)
    write(unit, *) 'lattice_dimension = ', Lat%dim
    write(unit, '(a, i15)') ' Lt    = ', ntime
    write(unit, '(a, f20.15)') ' dtau  = ', delt
    write(unit, '(a, f20.15)') ' Beta  = dtau*Lt = ', beta
    write(unit, '(a, f20.15)') ' lambda, el-phonon coupling = ', abs(ep_parameter)
    if(M > 0.0001) then
      write(unit, '(a, f20.15)') ' omega   = ', omega
      write(unit, '(a, f20.15)') ' M   = ', M
      write(unit, '(a, f20.15)') ' K   = ', D
      write(unit, '(a, f20.15)') ' U_eff = lam^2/omega^2   = ', U
    else
      write(unit, '(a)') 'Anti-aidabatic limit, M = 0.0'
    end if
    write(unit, '(a, f20.15)') ' disorder = ', disorder
    write(unit, '(a, f20.15)') ' filling  = ', filling
    write(unit, '(a, i15)') ' N_up = N_down  = ', nelec
    write(unit, '(a)') ''

    ! 写入观测值
    write(unit, '(a)') 'List of observables, averages and errors'
    write(unit, '(a)') repeat(' ', 1) // 'Description: '
    write(unit, '(a)') repeat(' ', 3) // 'e_PE = -U/2 * (n_up + n_down)^2 = -U/2*(n_up+n_down) - U*n_up*n_down'
    write(unit, '(a)') repeat(' ', 3) // 'e_ph_E = < g * X * (n_up + n_down)>'
    write(unit, '(a)') repeat(' ', 3) // 'For Hubbard model, e_TE = (e_KE + e_PE) '
    write(unit, '(a)') repeat(' ', 3) // 'For Holstein model, e_TE = (e_KE + e_ph_E) '
    write(unit, '(a)') repeat(' ', 3) // 'All observables are divided by Ns except N_tot'
    write(unit, '(a)') ''
    close(unit)
    
    ! create a .nml file to store the task_name of the info_sheet
    write(ci,'(a,1i4)') 'b',MPI_block
    ci = trim(output_addr)//'out_files/'//trim(ci)//'.nml'
    call remove_spaces(ci,len(ci))
    open(newunit=unit, file=ci, status='replace', action='write')
    write(unit, '(a)') '&file'
    write(unit, '(a)') ' filename = "'//trim(task_name)//'"'
    write(unit, '(a)') '/'
    close(unit)

  end subroutine

  subroutine remove_spaces(input_str,n)
    implicit none
    integer :: n
    character(len=*) :: input_str      ! 输入字符串
    character(len=n) :: output_str    ! 输出字符串（长度与输入相同）
    integer :: i, j

    output_str = ' '  ! 初始化为全空格
    j = 1             ! 指向输出字符串的当前位置

    do i = 1, len(input_str)
      if (input_str(i:i) /= ' ') then
        if (j > len(output_str)) exit  ! 防止越界
        output_str(j:j) = input_str(i:i)
        j = j + 1
      end if
    end do

    output_str = trim(output_str)  ! 去掉尾部空格
    input_str = output_str  ! 将结果赋值回输入字符串
  end subroutine remove_spaces

end MODULE input