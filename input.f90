
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
  integer :: MPI_nblock = 4!> = numprocs/MPI_one_block
! MC parameter
  integer :: warmup = 2, meas_interval = 1, meas_interval_tau = 20, meas_number = 0
  integer :: nbin_per_core = 20, nmeas_per_bin = 10! # of bins and size of one bin
! measurement parameter
  integer :: n_suit_corf = 8 , n_suit_ctau = 1 ! number of suited correlation functions
  logical :: file_exist
  character*30, allocatable :: obs_name(:),corf_name(:),ctau_name(:)
  integer,allocatable :: corf_length(:)
  real(8) :: max_pair = 0d0
  integer :: ntau = 80, one_bin_tau 

! General ARGUMENTS
  character(len=100) :: nml_file = 'input.nml', output_file = 'data/',output_format = 'f18.8'
  integer :: nflv = 1, ncopy = 1 ! number of flavours
  integer :: ntime, Ns ! ntime is number of time slices and Ns is # of sites
  integer :: nblock, ngroup = 5 !> nblock = ntime/ngroup + 2
  real(8) :: delt = 0.1, beta = 80d0,hop = 1d0
  integer :: print_loop = 2000
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

end MODULE input