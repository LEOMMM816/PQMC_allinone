module dqmc_measurements
  use pf_setting
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: MeasurementSystem

  integer, parameter :: dp = real64
  integer, parameter :: NAME_LEN = 64
  integer, parameter :: KIND_SCALAR = 1, KIND_FIELD = 2

  type :: MeasurementSystem
    ! ---- registry ----
    character(len=NAME_LEN), allocatable :: names(:)
    integer,                  allocatable :: kinds(:)      ! KIND_SCALAR / KIND_FIELD
    integer,                  allocatable :: widths(:)     ! 1 for scalar, field_width for field
    integer,                  allocatable :: off_lo(:), off_hi(:)
    integer :: nobs = 0
    integer :: total_width = 0
    integer :: field_width = 0
    integer,allocatable :: handles(:)
    ! ---- data buffer ----
    real(dp), allocatable :: data(:,:,:) ! [total_width, nbins, nmeas]
    integer :: nbins = 0, nmeas = 0
    integer :: cur_bin = 1, cur_meas = 0
    ! ---- status ----
    logical :: finalized = .false.

  contains
    ! ---- init / io ----
    procedure :: init_from_namelist   => ms_init_from_namelist   ! read observables, finalize, allocate data
    procedure :: finalize_layout      => ms_finalize_layout      ! compute offsets
    procedure :: init_storage         => ms_init_storage

    ! ---- run control ----
    procedure :: begin_measure            => ms_begin_measure
    ! ---- take measurement ----
    procedure :: take_measurement     => ms_take_measurement
    ! ---- record (handle only for speed) ----
    procedure :: get_handle           => ms_get_handle
    procedure :: get_range_by_handle  => ms_get_range_by_handle
    procedure :: record_scalar      => ms_record_scalar
    procedure :: record_field       => ms_record_field
    procedure :: record_field_entry => ms_record_field_entry

    ! ---- analysis ----
    procedure :: compute_all_bin_means => ms_compute_all_bin_means
  end type MeasurementSystem

  type(MeasurementSystem) :: Meas_sys
contains
!======================================================================!
!                             TAKE MEASUREMENTS                           !
!======================================================================!

  SUBROUTINE ms_take_measurement(this,forward, time)

    implicit none

    class(MeasurementSystem), intent(inout) :: this
    logical, intent(in) :: forward
    integer, intent(in) :: time
    integer :: head, bin , handle, N_cell
    integer ::   i,j,c1,c2,s1, s2, s1_next,s2_next,ind
    complex(kind=8) :: factor ! phase factor of the position (should be complex in general)
    real(8) :: obs_temp,  t_s1_c, t_s2_c
    type(cell_type), pointer :: pc1,pc2

      !! MEASUREMENT BEGINS
    head = this%cur_meas
    bin = this%cur_bin
    N_cell = Lat%N_cell

    do c1 = 1, N_cell
      ! preparation
      pc1 => p_cells(c1)
      s1 = pc1%sites(1)%id
      ! obs: phonon kinetic energy
      handle = this%get_handle('BF_KE')
      obs_temp = 0d0
      obs_temp = obs_temp + 0.5d0 * D * sum(boson_field(pc1%bf_list,time)**2) &
      & - biased_phonon * sum(boson_field(pc1%bf_list,time))
      obs_temp = obs_temp/N_cell
      call this%record_scalar(handle, obs_temp)

      !! obs: phonon kinetic energy
      handle = this%get_handle('BF_PE')
      obs_temp = 0d0
      obs_temp = obs_temp  - 0.5d0 * M * sum((boson_field(pc1%bf_list,time) &
      & - boson_field(pc1%bf_list,time+1))**2)/(delt)**2
      obs_temp = obs_temp/N_cell  + 1d0/(2d0*delt)
      call this%record_scalar(handle, obs_temp)

      !! obs: electron kinetic energy along x direction
      handle = this%get_handle('El_KEx')
      obs_temp = 0d0
      factor = (-1d0)**(1+(s1-1)/La)
      s1_next = s1 + 2
      if (mod(s1_next, La) == 1) s1_next = s1_next - La
      obs_temp = obs_temp + factor * (hop) * 2 * real(g_h(s1_next,s1) + g_h(s1,s1_next))
      obs_temp = obs_temp/N_cell
      call this%record_scalar(handle, obs_temp)

      !! obs : electron kinetic energy along y direction
      handle = this%get_handle('El_KEy')
      obs_temp = 0d0
      factor = (-1d0)**(1+(s1-1)/La)
      s1_next = s1 + La
      if(s1_next > Ns) s1_next = s1_next - Ns
      s1_next = s1_next + 1
      if (mod(s1_next, La) == 1) s1_next = s1_next - La
      obs_temp = obs_temp + factor * (-hop)  * 2 * real(g_h(s1_next,s1) + g_h(s1,s1_next))
      obs_temp = obs_temp/N_cell
      call this%record_scalar(handle, obs_temp)

      !! obs: electron potential energy
      handle = this%get_handle('El_PE')
      obs_temp = 0d0
      !obs_temp = obs_temp + 2 * real(ep_parameter * boson_field(s1,time) * (g_h(s1,s1)))
      !obs_temp = obs_temp -  real(ep_parameter**2 * g_h(s1,s1) * g_h(s1,s1)) - ep_parameter**2 * g_h(s1,s1)
      obs_temp = obs_temp  +  real(ep_parameter)**2 * g_h(s1,s1) * conjg(g_h(s1,s1))
      obs_temp = obs_temp/N_cell
      call this%record_scalar(handle, obs_temp)

      !! obs: electron-phonon energy
      handle = this%get_handle('ElPH_E')
      obs_temp = 0d0
      obs_temp = obs_temp + 2 * real(ep_parameter * boson_field(s1,time) * (g_h(s1,s1)))
      obs_temp = obs_temp/N_cell
      CALL this%record_scalar(handle, obs_temp)

      !! obs: electron total energy
      handle = this%get_handle('El_E')
      obs_temp = 0d0
      factor = (-1d0)**(1+(s1-1)/La)
      s1_next = s1 + 1
      if (mod(s1_next, La) == 1) s1_next = s1_next - La
      obs_temp = obs_temp + factor * (-hop) * 2 * real(g_h(s1_next,s1) + g_h(s1,s1_next))
      factor = (-1d0)**(1+(s1-1)/La)
      s1_next = s1 + La
      if(s1_next > Ns) s1_next = s1_next - Ns
      s1_next = s1_next + 1
      if (mod(s1_next, La) == 1) s1_next = s1_next - La
      obs_temp = obs_temp + factor * (-hop)  * 2 * real(g_h(s1_next,s1)+ g_h(s1,s1_next))
      if(M < 0.0001d0) then
        ! hubbard model
        obs_temp = obs_temp -  real(ep_parameter)**2 * g_h(s1,s1) * conjg(g_h(s1,s1))
      else
        !electron-phonon model
        obs_temp = obs_temp + 2 * real(ep_parameter * boson_field(s1,time) * (g_h(s1,s1)))
      end if
      obs_temp = obs_temp/N_cell
      CALL this%record_scalar(handle, obs_temp)

      !! obs : ph_X
      handle = this%get_handle('BF_X')
      obs_temp = 0d0
      obs_temp = obs_temp + (boson_field(s1,time))/Ns
      this%data(8,head, bin) = this%data(8,head, bin) + obs_temp
      obs_temp = obs_temp/N_cell
      call this%record_scalar(this%get_handle('BF_X'), obs_temp)

      !! obs: electron density
      handle = this%get_handle('El_den')
      obs_temp = 0d0
      obs_temp = obs_temp + 2 * real(g_h(s1,s1))
      this%data(9,head, bin) = this%data(9,head, bin) + obs_temp
      obs_temp = obs_temp/N_cell
      call this%record_scalar(handle, obs_temp)

      !! obs: double occupancy
      handle = this%get_handle('El_docc')
      obs_temp = 0d0
      obs_temp = obs_temp + (g_h(s1,s1) * conjg(g_h(s1,s1)))
      obs_temp = obs_temp/N_cell
      call this%record_scalar(handle, obs_temp)

      do c2 = 1, N_cell
      !! preparation for correlation functions
        pc2 => p_cells(c2)
        s1 = pc1%sites(1)%id
        s2 = pc2%sites(1)%id
        call get_relative_index(ind,pc1%dpos,pc2%dpos)

      !! corf1: A-A 1-particle
        handle = this%get_handle('corf_G1')
        obs_temp = 0d0
        obs_temp = obs_temp + g(s2,s1)
        obs_temp = obs_temp/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf2: A-A 2-particle
        handle = this%get_handle('corf_G2')
        obs_temp = 0d0
        obs_temp = obs_temp + g_h(s2,s1) * conjg(g_h(s2,s1))
        obs_temp = obs_temp/N_cell
        call this%record_field_entry(handle, ind, obs_temp)
      !! corf3: nu_nu
        handle = this%get_handle('corf_nunu')
        obs_temp = 0d0
        obs_temp = obs_temp + g_h(s2,s1) * g(s1,s2) + g_h(s1,s1) * (g_h(s2,s2))
        obs_temp = obs_temp/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf4: nu_nd
        handle = this%get_handle('corf_nund')
        obs_temp = 0d0
        obs_temp = obs_temp + g_h(s1,s1) * conjg(g_h(s2,s2))
        obs_temp = obs_temp/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf5: den-den
        handle = this%get_handle('corf_nn')
        obs_temp = 0d0
        obs_temp = obs_temp + (g_h(s1,s1) + conjg(g_h(s1,s1))) * (g_h(s2,s2) + conjg(g_h(s2,s2)))
        obs_temp = obs_temp + g_h(s2,s1) * g(s1,s2) + conjg(g_h(s2,s1) * g(s1,s2))
        obs_temp = obs_temp/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf6: X-X
        handle = this%get_handle('corf_XX')
        obs_temp = 0d0
        s1 =pc1%bf_list(1)
        s2 =pc2%bf_list(1)
        obs_temp = obs_temp + (boson_field(s1,time)) * (boson_field(s2,time))
        obs_temp = obs_temp/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      end do !for c2
    end do !for c1

      !! END MEASUREMENT

  END SUBROUTINE ms_take_measurement

!======================================================================!
!                             NAMELIST INIT                             !
!======================================================================!

  subroutine init_meas()
    implicit none
    Meas_sys%nbins = nbin_per_core
    Meas_sys%nmeas = n_meas_per_bin
    Meas_sys%cur_bin = 1
    Meas_sys%cur_meas = 1
    call Meas_sys%init_from_namelist('measurement.nml')
  end subroutine init_meas

  subroutine ms_init_from_namelist(this, filename)
    class(MeasurementSystem), intent(inout) :: this
    character(*),            intent(in)    :: filename
    integer :: nbins, nmeas_per_bin
    integer :: lun, ios
    character(len=NAME_LEN) :: name, mskind
    integer :: cap

    namelist /observable/ name, mskind

    if (nbins<=0 .or. nmeas_per_bin<=0) error stop "init_from_namelist: nbins/nmeas must be >0"
    this%field_width = product(Lat%dlength)

    ! 初始容量（动态增长，避免频繁重分配）
    cap = 16
    allocate(this%names(cap), this%kinds(cap), this%widths(cap), this%off_lo(cap), this%off_hi(cap))
    this%nobs = 0

    open(newunit=lun, file=filename, status='old', action='read', iostat=ios)
    if (ios/=0) error stop "init_from_namelist: cannot open file: "//trim(filename)

    do
      name = ''; mskind = ''
      read(lun, nml=observable, iostat=ios)
      if (ios /= 0) exit  ! 到末尾或无更多组
      call tolower_inplace(mskind)
      call trim_inplace(name)
      if (len_trim(name)==0) error stop "observable group missing name"
      select case (mskind)
      case ('scalar','s','sc','1')
        call push_obs(this, name, KIND_SCALAR)
      case ('field','f','array','a')
        if (this%field_width<=0) error stop "observable mskind=field but field_width not provided"
        call push_obs(this, name, KIND_FIELD)
      case default
        error stop "unknown observable mskind: "//trim(mskind)//" (name="//trim(name)//")"
      end select
    end do
    close(lun)

    if (this%nobs==0) error stop "no observable groups found in namelist"

    call this%finalize_layout()
    call this%init_storage(nbins, nmeas_per_bin)
  end subroutine ms_init_from_namelist

  subroutine push_obs(this, name, k)
    class(MeasurementSystem), intent(inout) :: this
    character(*),            intent(in)    :: name
    integer,                 intent(in)    :: k
    integer :: i

    ! 去重（线性查找，仅初始化时调用一次，无性能瓶颈）
    do i=1,this%nobs
      if (trim(this%names(i)) == trim(name)) error stop "duplicate observable name: "//trim(name)
    end do

    ! 扩容
    if (this%nobs == size(this%names)) then
      call grow_arrays(this)
    end if

    this%nobs = this%nobs + 1
    this%names(this%nobs) = adjustl(name)
    this%kinds(this%nobs) = k
    select case (k)
    case (KIND_SCALAR); this%widths(this%nobs) = 1
    case (KIND_FIELD); this%widths(this%nobs) = this%field_width
    end select
  end subroutine push_obs

  subroutine grow_arrays(this)
    class(MeasurementSystem), intent(inout) :: this
    integer :: newcap
    character(len=NAME_LEN), allocatable :: n2(:)
    integer, allocatable :: i2(:)

    newcap = max(16, int(1.5d0*size(this%names)) + 16)

    allocate(n2(newcap)); n2(1:this%nobs) = this%names(1:this%nobs); call move_alloc(n2, this%names)
    allocate(i2(newcap)); i2(1:this%nobs) = this%kinds(1:this%nobs); call move_alloc(i2, this%kinds)
    allocate(i2(newcap)); i2(1:this%nobs) = this%widths(1:this%nobs); call move_alloc(i2, this%widths)
    allocate(i2(newcap)); i2(1:this%nobs) = this%off_lo(1:this%nobs); call move_alloc(i2, this%off_lo)
    allocate(i2(newcap)); i2(1:this%nobs) = this%off_hi(1:this%nobs); call move_alloc(i2, this%off_hi)
  end subroutine grow_arrays

  subroutine tolower_inplace(s)
    character(len=*), intent(inout) :: s
    integer :: i, ia
    do i=1,len_trim(s)
      ia = iachar(s(i:i))
      if (ia>=iachar('A') .and. ia<=iachar('Z')) s(i:i) = achar(ia+32)
    end do
  end subroutine tolower_inplace

  subroutine trim_inplace(s)
    character(len=*), intent(inout) :: s
    s = trim(adjustl(s))
  end subroutine trim_inplace

!======================================================================!
!                            LAYOUT & STORAGE                           !
!======================================================================!

  subroutine ms_finalize_layout(this)
    class(MeasurementSystem), intent(inout) :: this
    integer :: i, cursor
    if (this%nobs<=0) error stop "finalize_layout: no observables"
    cursor = 1
    do i=1,this%nobs
      this%off_lo(i) = cursor
      this%off_hi(i) = cursor + this%widths(i) - 1
      cursor = this%off_hi(i) + 1
    end do
    this%total_width = cursor - 1
    this%finalized = .true.
  end subroutine ms_finalize_layout

  subroutine ms_init_storage(this, nbins, nmeas_per_bin)
    class(MeasurementSystem), intent(inout) :: this
    integer,                 intent(in)    :: nbins, nmeas_per_bin
    if (.not.this%finalized) error stop "init_storage: call finalize_layout first"
    if (nbins<=0 .or. nmeas_per_bin<=0) error stop "init_storage: nbins/nmeas must be >0"
    this%nbins = nbins
    this%nmeas = nmeas_per_bin
    allocate(this%data(this%total_width, this%nmeas,this%nbins))
    this%data = 0.0_dp
    this%cur_bin  = 1
    this%cur_meas = 0
  end subroutine ms_init_storage

!======================================================================!
!                         RUN CONTROL & RECORD                          !
!======================================================================!
  subroutine ms_begin_measure(this,forward, time)
    class(MeasurementSystem), intent(inout) :: this
    logical, intent(in) :: forward
    integer, intent(in) :: time
    !2nd trotter
    if(second_trotter) then

    end if
    call take_measurement(this, time)
    ! next loop
    if ((forward .and. time == ntime/2 + meas_number) .or. ((.not. forward) .and. time == ntime/2 - meas_number)) then
      ! average over all sites and all measurements in one loop
      this%data(:,this%cur_meas, this%cur_bin) = this%data(:,this%cur_meas, this%cur_bin)/(2*meas_number + 1)
      this%cur_meas = this%cur_meas + 1
      if (this%cur_meas > this%nmeas) then
        ! next bin
        this%cur_meas = 1
        this%cur_bin = this%cur_bin + 1
        if (this%cur_bin > this%nbins) error stop "begin_measure: no more bins left"
      end if

    end if

    ! next loop

    !2nd trotter
    if(second_trotter) then

    end if
  end subroutine ms_begin_measure

  integer function ms_get_handle(this, name) result(h)
    class(MeasurementSystem), intent(in) :: this
    character(*),            intent(in) :: name
    integer :: i
    do i=1,this%nobs
      if (trim(this%names(i)) == trim(name)) then
        h = i; return
      end if
    end do
    h = 0
    if (h==0) error stop "get_handle: name not found: "//trim(name)
  end function ms_get_handle

  subroutine ms_get_range_by_handle(this, handle, lo, hi, width, mskind)
    class(MeasurementSystem), intent(in)  :: this
    integer,                 intent(in)  :: handle
    integer,                 intent(out) :: lo, hi, width, mskind
    if (handle<1 .or. handle>this%nobs) error stop "get_range_by_handle: handle out of range"
    lo = this%off_lo(handle); hi = this%off_hi(handle)
    width = this%widths(handle); mskind  = this%kinds(handle)
  end subroutine ms_get_range_by_handle

  subroutine ms_record_scalar(this, handle, value, bin, meas)
    class(MeasurementSystem), intent(inout) :: this
    integer,                 intent(in)    :: handle
    real(dp),                intent(in)    :: value
    integer,      intent(in), optional    :: bin, meas
    integer :: lo, hi, width, mskind, b, m
    call this%get_range_by_handle(handle, lo, hi, width, mskind)
    if (mskind /= KIND_SCALAR) error stop "record_scalar: not a scalar handle"
    b = merge(bin, this%cur_bin, present(bin))
    m = merge(meas, this%cur_meas, present(meas))
    if (b<1 .or. b>this%nbins) error stop "record_scalar: bin out of range"
    if (m<1 .or. m>this%nmeas) error stop "record_scalar: meas out of range"
    this%data(lo, m, b) = this%data(lo, m, b) + value

  end subroutine ms_record_scalar

  subroutine ms_record_field(this, handle, arr, bin, meas)
    class(MeasurementSystem), intent(inout) :: this
    integer,                 intent(in)    :: handle
    real(dp),                intent(in)    :: arr(:)   ! flattened field (length = field_width)
    integer,      intent(in), optional    :: bin, meas
    integer :: lo, hi, width, mskind, b, m
    call this%get_range_by_handle(handle, lo, hi, width, mskind)
    if (mskind /= KIND_FIELD) error stop "record_field: not a field handle"
    if (size(arr) /= width) error stop "record_field: size(arr)!=field width"
    b = merge(bin, this%cur_bin, present(bin))
    m = merge(meas, this%cur_meas, present(meas))
    if (b<1 .or. b>this%nbins) error stop "record_field: bin out of range"
    if (m<1 .or. m>this%nmeas) error stop "record_field: meas out of range"
    this%data(lo:hi, m,b) = this%data(lo:hi, m,b) + arr

  end subroutine ms_record_field

  subroutine ms_record_field_entry(this, handle, idx, value, bin, meas)
    class(MeasurementSystem), intent(inout) :: this
    integer,                 intent(in)    :: handle
    integer,                 intent(in)    :: idx       ! index in the flattened field (1..field_width)
    real(dp),                intent(in)    :: value
    integer,      intent(in), optional    :: bin, meas
    integer :: lo, hi, width, mskind, b, m
    call this%get_range_by_handle(handle, lo, hi, width, mskind)
    if (mskind /= KIND_FIELD) error stop "record_field_entry: not a field handle"
    if (idx<1 .or. idx>width) error stop "record_field_entry: idx out of range"
    b = merge(bin, this%cur_bin, present(bin))
    m = merge(meas, this%cur_meas, present(meas))
    if (b<1 .or. b>this%nbins) error stop "record_field_entry: bin out of range"
    if (m<1 .or. m>this%nmeas) error stop "record_field_entry: meas out of range"
    this%data(lo+idx-1, m, b) = this%data(lo+idx-1, m, b) + value
  end subroutine ms_record_field_entry
!======================================================================!
!                               ANALYSIS                                !
!======================================================================!
  subroutine ms_data_analysis(this)
    class(MeasurementSystem), intent(inout) :: this
    ! 这里可以添加更多分析功能
    
  end subroutine ms_data_analysis
  subroutine ms_compute_all_bin_means(this, means)
    ! 计算所有观测量在每个 bin 的均值（沿第 3 维求平均）
    ! 返回: means(total_width, nbins)
    class(MeasurementSystem), intent(in)  :: this
    real(dp), allocatable,    intent(out) :: means(:,:)
    integer :: b, m

    if (.not.allocated(this%data)) error stop "compute_all_bin_means: storage not initialized"

    allocate(means(this%total_width, this%nbins))
    means = 0.0_dp
    do b = 1, this%nbins
      do m = 1, this%nmeas
        means(:, b) = means(:, b) + this%data(:, m, b)
      end do
      means(:, b) = means(:, b) / real(this%nmeas, dp)
    end do
  end subroutine ms_compute_all_bin_means

!======================================================================!
!                        USER MEASUREMENT PLACEHOLDER                   !
!======================================================================!
!  你可以在主程序里这样使用：
!    call M%start_bin(b)
!    do s=1,nmeas_per_bin
!      call M%begin_measurement()
!      ! >>> YOUR MEASUREMENT CODE HERE <<<
!      ! 例如：
!      !   call M%record_scalar(h_energy, energy_value)
!      !   call M%record_field (h_green,  gf_flattened(:))
!    end do
!
!  建议在初始化后一次性用 get_handle 拿到所有需要的 handle：
!    h_energy = M%get_handle("energy")
!    h_green  = M%get_handle("green_fn")
!======================================================================!

end module dqmc_measurements
