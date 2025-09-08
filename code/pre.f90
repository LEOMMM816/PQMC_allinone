# 1 "input.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "input.f90"

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
  logical :: global_update = .true.
  logical :: global_update_shift = .false.
  logical :: global_update_flip = .true.
  logical :: global_update_exchange = .false.
  logical :: record_middata = .false.
  logical,target :: debug = .false.
  logical :: meas_corf = .true.
  logical :: fixedseeds = .false.
  logical :: checkerboard = .true.
  logical :: second_trotter = .false.

! constant
  real(8),parameter :: e = 2.7182818284590452353602874713526624d0
  real(8),parameter :: pi = 3.1415926535897932384626433832795029d0
  character*30, parameter :: file_name = 'data/'
! mpi parameter
  integer :: ierr, myid = 0, numprocs = 1, MPI_block = 1
  integer :: MPI_one_block = 8 !> number of cores in one mpi_block, = 1 if not mpi
  integer :: MPI_nblock = 4!> = numprocs/MPI_one_block
! MC parameter
  integer :: warmup = 20000, meas_interval = 1, meas_interval_tau = 20, meas_number = 0
  integer :: nbin_per_core = 20, one_bin = 2000! # of bins and size of one bin
! measurement parameter
  integer :: n_obs = 11 , n_corf = 33, n_ctau = 5 ! number of observables
  integer :: n_suit_corf = 8 , n_suit_ctau = 1 ! number of suited correlation functions
  integer :: data_per_line = 100
  logical :: file_exist
  character*30, allocatable :: obs_name(:),corf_name(:),ctau_name(:)
  integer :: corf_length
  real(8) :: max_pair = 0d0
  integer :: ntau = 80, one_bin_tau

! General ARGUMENTS
  integer :: nflv = 1, ncopy = 1 ! number of flavours
  integer :: La =40, Lb = 2! one row has La sites and one column has Lb sites (La >= Lb >= 1,La > 2,2*2->1*4)
  integer :: ntime, Ns ! ntime is number of time slices and Ns is # of sites
  integer :: nblock, ngroup = 5 !> nblock = ntime/ngroup + 2
  real(8) :: delt = 0.1, beta = 80d0
  integer :: print_loop = 5000
  integer :: lattice_dimension
  real(8) :: prob !> accept probobility
  real(8) :: err_fast = 0d0, err_fast_max = 0.000001d0
  integer, allocatable :: field_type(:) ! field number
  real(8):: TBC_phase_number(2) = (/0d0,0d0/)
  real(8) :: filling
  real(8) :: biased_phonon

! HS parameter
  integer :: isingfield_number = 1 ! type of HS fields
  real(kind=8) :: U = 0d0 ! Hubbard U and steps delt * U < 0.5
  real(kind=8) :: miu  ! chemical potential
  integer, allocatable :: isingfieldmax(:) ! (isingfield_number)
  integer, allocatable :: isingfield(:, :, :) ! discrete ising field space*time*isingfield_number
! proj
  real(8) :: disorder = 0d0
  real(kind=8) :: eta ! et a = exp(-delt*(miu - U/2))
  complex(kind=8) :: hop = 1d0 ! hopping term
  integer,parameter :: cb_site_num = 4, cb_bond_num =4
  integer :: n_cb_pla !> number of plaquettes in checkerboard
  integer :: cb_field_number  !> number of decomposed checkerboard field
  real(8) :: filling_factor
  integer :: nelec
  real(kind=8) R_nflv ! relative boltzmann weight for different particles
  complex(8), ALLOCATABLE :: slater(:, :), slater_Q(:, :), slater_D(:)
  complex(kind=8), allocatable :: K_slater(:, :)! kinetic energy ns*ns

! phonon parameter
  integer :: n_phonon_field = 4 !> number of phonon fields
  integer :: n_bf = 1 ! number of boson fields
  integer :: phonon_field_number ! number of decomposed phonon field
  real(8) :: D = 1d0,M = 1d0 !>stiffness constant and Mass
  real(8) :: char_length !> characteristic length of phonon field
  complex(8) :: ep_parameter!>electron-phonon coupling
  real(8) :: omega!> phonon frequency
  real(8) :: end_field = 0d0
  real(8) :: jump_distance
  real(8) :: max_displacement !> = hop / ep_parameter
  integer :: n_local_update = 1 !> times of proposal  of new ph field at each location
  real(8), allocatable :: boson_field(:,:) !(n_bf*N_cell,time)


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
    integer,allocatable :: bf_list(:) !> index of the coupled boson field
    complex,allocatable :: expKV(:,:,:,:) !> exp(-delt* K or V) for each pf, (ppf%dim,ppf%dim,n_plaquette,ntime)
    complex,allocatable :: expKV_inv(:,:,:,:) !> exp(delt* K) for each pf, (ppf%dim,ppf%dim,n_plaquette,ntime)
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
    logical :: updateable
    integer :: dim
    integer :: n_plaquette
    complex(8),allocatable :: Vmatrix(:,:)
    complex(8),allocatable :: Kmatrix(:,:)
    real(8), allocatable :: pla_tsl_vec(:,:) ! posotions of phonon field, translational vector,(lat%dim,lat%dim)
    real(8), allocatable :: pla_tsl_offset_vec(:) ! offset of phonon field, translational vector (lat%dim)
    real(8), allocatable :: pla_int_vec(:,:) ! positions of sites in each plaquette, relative to the first site,(lat%dim,ppf%dim)
    type(pf_data_type),pointer :: p_pfdata !> phonon field data
    logical,pointer :: debug
  end type pf_type
  
 
  type:: cell_type
    integer :: id! unique cell id
    integer :: n_site,n_bf ! number of sites and boson fields in this cell
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


contains

  SUBROUTINE init(block)
    implicit none
    integer, intent(in) :: block
    integer ::  n ,i,j

    ! parameters definition

    one_bin_tau = one_bin/(meas_interval_tau/meas_interval)
    !print*,'one_bin_tau:',one_bin_tau
    call init_ep_parameter(block)

    jump_distance = (1d0) * (0.1d0/delt)

    ntime = nint(beta/delt)
    delt = beta/ntime
    Ns = Lat%Ns
    nblock = ntime/ngroup + 2
    D = 1d0 * omega**2
    U = real(ep_parameter)**2 /D
    ! phonon :set the bond and plaquette
    corf_length = Lat%dlength(1)
    phonon_field_number = 1
    
    filling_factor = filling
    nelec = nint(Ns*filling_factor)
    
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
  end subroutine

  subroutine init_ep_parameter(mpi_block)

    integer, intent(in) :: mpi_block
    ep_parameter = 4d0
    ep_parameter = complex(0d0,sqrt(real(ep_parameter)))
    omega = sqrt(2d0)**(mpi_block)
    filling = 10d0/40d0

    return
  end subroutine init_ep_parameter

end MODULE input
