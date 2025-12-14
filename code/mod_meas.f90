module dqmc_measurements
  use pf_setting
  implicit none
  public :: MeasurementSystem

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
    ! ---- data buffer ----
    complex(dp), allocatable :: data(:,:,:) ! [total_width, nbins, nmeas]
    integer :: nbins = 0, nmeas = 0
    integer :: cur_bin = 1, cur_meas = 0
    ! ---- status ----
    logical :: finalized = .false.
    ! ---- statistics ----
    real(dp), allocatable :: mean(:,:)  ! mean for each observable (total_width)

  contains
    ! ---- init / io ----
    procedure :: init_meas            => ms_init_meas               ! initialize measurement system
    procedure :: init_from_namelist   => ms_init_from_namelist   ! read observables, finalize, allocate data
    procedure :: finalize_layout      => ms_finalize_layout      ! compute offsets
    procedure :: init_storage         => ms_init_storage

    ! ---- run control ----
    procedure :: begin_measure        => ms_begin_measure
    ! ---- take measurement ----
    procedure :: take_measurement     => ms_take_measurement
    procedure :: k_space_obs       => ms_k_space_obs
    ! ---- record (handle only for speed) ----
    procedure :: get_handle           => ms_get_handle
    procedure :: get_range_by_handle  => ms_get_range_by_handle
    procedure :: record_scalar      => ms_record_scalar
    procedure :: record_field       => ms_record_field
    procedure :: record_field_entry => ms_record_field_entry

    ! ---- analysis ----
    procedure :: data_analyse         => ms_data_analyse
    procedure :: compute_all_bin_means => ms_compute_all_bin_means
    ! ---- output ----
    procedure :: output_array_data => ms_output_array_data
    procedure :: output_reading_guide => ms_output_reading_guide
  end type MeasurementSystem

  type(MeasurementSystem) :: Meas_sys
contains
!======================================================================!
!                             TAKE MEASUREMENTS                           !
!======================================================================!

!=========================tool functions for measurements========================!
  subroutine spin_current_matrix(pc1,pc2,pc3,pc4,mat,n_g,hatree_only)
    ! this function computes the two-body matrix element of the spin current operator
    ! spin current correation between i-j bond and k-l bond
    ! [ci1* cj2 , ci2* cj1, cj1* ci2, cj2* ci1] [i,j -> k,l]
    ! pc1,pc2 are the left operators, pc3,pc4 are the right operators, 1,2 are spin indices
    ! n_g is the number of one-body operators in each current operator, for spin_x(y) current n_g = 4
    ! pc1%sites(1) = i1, pc1%sites(2) = i2, pc2%sites(1) = j1, pc2%sites(2) = j2, etc.
    ! to retrieve spin_x(y) current, one need to contract the matrix with vecters [-i,-i,i,i] ([-1,1,1,-1])
    implicit none
    integer, intent(in) :: n_g ! # of one-body operators in each current operator
    type(cell_type), pointer,intent(in) :: pc1,pc2,pc3,pc4
    complex(8), intent(inout) :: mat(n_g,n_g)
    integer :: ls(n_g,2), rs(n_g,2),i,j
    logical, intent(in) :: hatree_only
    mat = 0.0d0
    ls(1,:) = [pc1%sites(1),pc2%sites(2)]
    ls(2,:) = [pc1%sites(2),pc2%sites(1)]
    ls(3,:) = [pc2%sites(1),pc1%sites(2)]
    ls(4,:) = [pc2%sites(2),pc1%sites(1)]
    rs(1,:) = [pc3%sites(1),pc4%sites(2)]
    rs(2,:) = [pc3%sites(2),pc4%sites(1)]
    rs(3,:) = [pc4%sites(1),pc3%sites(2)]
    rs(4,:) = [pc4%sites(2),pc3%sites(1)]
    call cal_wick_contraction(mat, n_g, ls, rs, hatree_only)

  end subroutine spin_current_matrix

  subroutine spin_momentum_matrix(pc1,pc2,mat,n_g,hatree_only)
    ! this function computes the two-body matrix element of the spin momentum operator
    ! spin momentum correation between unit cell pc1 and pc2
    ! [ci1* ci1 , ci1* ci2, ci2* ci1, ci2* ci2] [i,j -> k,l]
    implicit none
    integer, intent(in) :: n_g != 4 in this case,# of one-body operators in each momentum operator
    type(cell_type), pointer,intent(in) :: pc1,pc2
    complex(8), intent(inout) :: mat(n_g,n_g)
    integer :: ls(n_g,2), rs(n_g,2),i,j
    logical, intent(in) :: hatree_only
    mat = 0.0d0
    ls(1,:) = [pc1%sites(1),pc1%sites(1)] ! 1_up, 1_up
    ls(2,:) = [pc1%sites(1),pc1%sites(2)] ! 1_up, 1_dn
    ls(3,:) = [pc1%sites(2),pc1%sites(1)] ! 1_dn, 1_up
    ls(4,:) = [pc1%sites(2),pc1%sites(2)] ! 1_dn, 1_dn
    rs(1,:) = [pc2%sites(1),pc2%sites(1)] ! 2_up, 2_up
    rs(2,:) = [pc2%sites(1),pc2%sites(2)] ! 2_up, 2_dn
    rs(3,:) = [pc2%sites(2),pc2%sites(1)] ! 2_dn, 2_up
    rs(4,:) = [pc2%sites(2),pc2%sites(2)] ! 2_dn, 2_dn

    call cal_wick_contraction(mat, n_g, ls, rs, hatree_only)
  end subroutine spin_momentum_matrix

  subroutine kinetic_matrix(pc1,pc2,pc3,pc4,mat,n_g,hatree_only)
    implicit none
    integer, intent(in) :: n_g ! # of one-body operators in each kinetic operator
    type(cell_type), pointer,intent(in) :: pc1,pc2,pc3,pc4
    complex(8), intent(inout) :: mat(n_g,n_g)
    integer :: ls(n_g,2), rs(n_g,2),i,j
    logical, intent(in) :: hatree_only
    mat = 0.0d0
    ls(1,:) = [pc1%sites(1),pc2%sites(1)]
    ls(2,:) = [pc2%sites(1),pc1%sites(1)]
    ls(3,:) = [pc1%sites(2),pc2%sites(2)]
    ls(4,:) = [pc2%sites(2),pc1%sites(2)]
    rs(1,:) = [pc3%sites(1),pc4%sites(1)]
    rs(2,:) = [pc4%sites(1),pc3%sites(1)]
    rs(3,:) = [pc3%sites(2),pc4%sites(2)]
    rs(4,:) = [pc4%sites(2),pc3%sites(2)]
    call cal_wick_contraction(mat, n_g, ls, rs, hatree_only)
  end subroutine kinetic_matrix

  subroutine physical_spin_current_matrix(pc1,pc2,pc3,pc4,mat,hatree_only)
    implicit none
    type(cell_type), pointer,intent(in) :: pc1,pc2,pc3,pc4
    complex(8), intent(inout) :: mat(8,8)
    integer :: ls(8,2), rs(8,2),i,j
    logical, intent(in) :: hatree_only
    mat = 0.0d0
    ! the first four are for spin-x -y current operator, the last four are for kinetic operator
    ls(1,:) = [pc1%sites(1),pc2%sites(2)]
    ls(2,:) = [pc1%sites(2),pc2%sites(1)]
    ls(3,:) = [pc2%sites(1),pc1%sites(2)]
    ls(4,:) = [pc2%sites(2),pc1%sites(1)]
    rs(1,:) = [pc3%sites(1),pc4%sites(2)]
    rs(2,:) = [pc3%sites(2),pc4%sites(1)]
    rs(3,:) = [pc4%sites(1),pc3%sites(2)]
    rs(4,:) = [pc4%sites(2),pc3%sites(1)]
    ls(5,:) = [pc1%sites(1),pc2%sites(1)]
    ls(6,:) = [pc2%sites(1),pc1%sites(1)]
    ls(7,:) = [pc1%sites(2),pc2%sites(2)]
    ls(8,:) = [pc2%sites(2),pc1%sites(2)]
    rs(5,:) = [pc3%sites(1),pc4%sites(1)]
    rs(6,:) = [pc4%sites(1),pc3%sites(1)]
    rs(7,:) = [pc3%sites(2),pc4%sites(2)]
    rs(8,:) = [pc4%sites(2),pc3%sites(2)]
    call cal_wick_contraction(mat, 8, ls, rs, hatree_only)

  end subroutine physical_spin_current_matrix

  subroutine cal_wick_contraction(mat, n_g, left_sites, right_sites, hatree_only)
    implicit none
    integer, intent(in) :: n_g
    complex(8), intent(inout) :: mat(n_g,n_g)
    integer, intent(in) :: left_sites(n_g,2), right_sites(n_g,2)
    logical, intent(in) :: hatree_only
    integer :: i,j
    mat = 0.0d0
    do i = 1,n_g
      do j = 1,n_g
        ! hatree part + fock part
        ! <ci1* cj2 , ci2* cj1> = <ci1* cj2> <ci2* cj1> + <ci1* cj1> <cj2* ci2>
        mat(i,j) = mat(i,j) + (g_h(left_sites(i,2),left_sites(i,1)) * g_h(right_sites(j,2),right_sites(j,1)))   ! hatree part
        if(.not. hatree_only) then
          mat(i,j) = mat(i,j) + g_h(right_sites(j,2),left_sites(i,1)) * g(left_sites(i,2),right_sites(j,1)) ! fork part
        end if
      end do
    end do
  end subroutine cal_wick_contraction
!=========================main measurement subroutine========================!
  SUBROUTINE ms_take_measurement(this, time)
    implicit none
    class(MeasurementSystem), intent(inout) :: this
    integer, intent(in) :: time
    integer :: head, bin , handle, N_cell,La,Lb,ind
    integer ::   i,j,c1,c2,s1_u, s1_d,bf1_x,bf1_y, s2_u, s2_d,bf2_x,bf2_y
    complex(kind=8) :: factor ! phase factor of the position (should be complex in general)
    complex(dp) :: obs_temp
    type(cell_type), pointer :: pc1,pc1_x,pc1_y, pc2, pc2_x,pc2_y
    complex(8), allocatable :: vec_l(:), vec_r(:)
    complex(8) :: spinJxy_mat(4,4),spinM_mat(4,4)
    complex(8) :: phy_spinJxy_mat(8,8)
    La = Lat%dlength(1)
    Lb = Lat%dlength(2)
    N_cell = Lat%N_cell
    head = this%cur_meas
    bin = this%cur_bin
    ! prepare G for measurement
    do j = 1, ns
          do i = 1, ns
            g(i, j) = -g_h(i, j)
          end do
          g(j, j) = g(j, j) + 1.0d0
    end do
      !! MEASUREMENT BEGINS
    !print*,"Taking measurement at time slice ", time, " bin ", bin, " meas ", head
    do c1 = 1, N_cell
    !! preparation
      pc1 => p_cells(c1)
      call get_uc_index_from_dpos(pc1%dpos + [1,0],ind)
      pc1_x => p_cells(ind)
      call get_uc_index_from_dpos(pc1%dpos + [0,1],ind)
      pc1_y => p_cells(ind)
      s1_u = pc1%sites(1)
      s1_d = pc1%sites(2)
      bf1_x = pc1%bf_list(1)
      bf1_y = pc1%bf_list(2)

    !! obs: phonon kinetic energy
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
      obs_temp = obs_temp + 1d0/(2d0*delt) * 2
      obs_temp = obs_temp/N_cell
      call this%record_scalar(handle, obs_temp)
    !! obs : ph_X
      handle = this%get_handle('BF_X')
      obs_temp = 0d0
      obs_temp = obs_temp + sum(boson_field(pc1%bf_list,time))/N_cell
      call this%record_scalar(this%get_handle('BF_X'), obs_temp)
    !! corfs are measured between two unit cells, pc2 is introduced
      do c2 = 1, N_cell
      !! preparation for correlation functions
        pc2 => p_cells(c2)
        call get_uc_index_from_dpos(pc2%dpos + [1,0],ind)
        pc2_x => p_cells(ind)
        call get_uc_index_from_dpos(pc2%dpos + [0,1],ind)
        pc2_y => p_cells(ind)
        s2_u = pc2%sites(1)
        s2_d = pc2%sites(2)
        bf2_x = pc2%bf_list(1)
        bf2_y = pc2%bf_list(2)
        call get_relative_index(ind,pc1%dpos,pc2%dpos)
        !print*,"Measuring corfs between cell ", c1, " and cell ", c2, " with relative index ", ind
        !print*,"  pc1 dpos ", pc1%dpos, " pc2 dpos ", pc2%dpos
        !print*,"  pc1_x's ind  ", pc1_x%id, " pc1_y's ind ", pc1_y%id
        !print*,"  pc2_x's ind  ", pc2_x%id, " pc2_y's ind ", pc2_y%id
        ! ind is the relative index of pc2 to pc1, used to index the correlation functions as an array data
        if(allocated(vec_l)) deallocate(vec_l, vec_r)
        allocate(vec_l(4), vec_r(4))
      !! corf: bf1-bf1
        handle = this%get_handle('BFx_BFx')
        obs_temp = 0d0
        obs_temp = obs_temp + (boson_field(bf1_x,time)) * (boson_field(bf2_x,time))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)
      !! corf: bf1-bf2
        handle = this%get_handle('BFx_BFy')
        obs_temp = 0d0
        obs_temp = obs_temp + (boson_field(bf1_x,time)) * (boson_field(bf2_y,time))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)
      !! corf: bf2-bf2
        handle = this%get_handle('BFy_BFy')
        obs_temp = 0d0
        obs_temp = obs_temp + (boson_field(bf1_y,time)) * (boson_field(bf2_y,time))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)
      !! corf: bf12-bf12
        handle = this%get_handle('BFxy_BFxy')
        obs_temp = 0d0
        obs_temp = obs_temp + (boson_field(bf1_x,time) * boson_field(bf1_y,time) & 
        & * boson_field(bf2_x,time) * boson_field(bf2_y,time))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)
      !! corf1: x-bond's spin-x current correlation
        call spin_current_matrix(pc1,pc1_x,pc2,pc2_x,spinJxy_mat,4,hatree_only = .false.) ! x-bond's spinJ mat
        handle = this%get_handle('Jxx_Jxx')
        obs_temp = 0d0
        vec_l = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf2: x-bond's spin-y current correlation
        handle = this%get_handle('Jxy_Jxy')
        obs_temp = 0d0
        vec_l = -hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf3: y-bond's spin-x current correlation
        call spin_current_matrix(pc1,pc1_y,pc2,pc2_y,spinJxy_mat,4,hatree_only = .false.) ! y-bond's spinJ mat
        handle = this%get_handle('Jyx_Jyx')
        obs_temp = 0d0
        vec_l = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf4: y-bond's spin-y current correlation
        handle = this%get_handle('Jyy_Jyy')
        obs_temp = 0d0
        vec_l = -hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_l))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf： x-bond's spin-x current and y-bond's spin-x current correlation
        call spin_current_matrix(pc1,pc1_x,pc2,pc2_y,spinJxy_mat,4,hatree_only = .false.) ! x-bond's spinJ mat
        handle = this%get_handle('Jxx_Jyx')
        obs_temp = 0d0
        vec_l = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        vec_r = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell

        call this%record_field_entry(handle, ind, obs_temp)
      !! corf： x-bond's spin-y current and y-bond's spin-y current correlation
        handle = this%get_handle('Jxy_Jyy')
        !! corf1: x-bond's spin-x current correlation
        call spin_current_matrix(pc1,pc1_x,pc2,pc2_x,spinJxy_mat,4,hatree_only = .true.) ! x-bond's spinJ mat
        handle = this%get_handle('Jxxh_Jxxh')
        obs_temp = 0d0
        vec_l = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf2: x-bond's spin-y current correlation
        handle = this%get_handle('Jxyh_Jxyh')
        obs_temp = 0d0
        vec_l = -hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf3: y-bond's spin-x current correlation
        call spin_current_matrix(pc1,pc1_y,pc2,pc2_y,spinJxy_mat,4,hatree_only = .true.) ! y-bond's spinJ mat
        handle = this%get_handle('Jyxh_Jyxh')
        obs_temp = 0d0
        vec_l = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf4: y-bond's spin-y current correlation
        handle = this%get_handle('Jyyh_Jyyh')
        obs_temp = 0d0
        vec_l = -hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_l))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf： x-bond's spin-x current and y-bond's spin-x current correlation
        call spin_current_matrix(pc1,pc1_x,pc2,pc2_y,spinJxy_mat,4,hatree_only = .true.) ! x-bond's spinJ mat
        handle = this%get_handle('Jxxh_Jyxh')
        obs_temp = 0d0
        vec_l = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        vec_r = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell

        call this%record_field_entry(handle, ind, obs_temp)
      !! corf： x-bond's spin-y current and y-bond's spin-y current correlation
        handle = this%get_handle('Jxyh_Jyyh')
        obs_temp = 0d0
        vec_l = -hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        vec_r = -hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)
      !! corf: x-bond's spin-x current and y-bond's spin-y current correlation
        call spin_current_matrix(pc1,pc1_x,pc2,pc2_y,spinJxy_mat,4,hatree_only = .true.) ! x-bond and y-bond
        handle = this%get_handle('Jxxh_Jyyh')
        obs_temp = 0d0
        vec_l = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        vec_r = -hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)
      !! corf: x-bond's spin-y current and y-bond's spin-x current correlation
        handle = this%get_handle('Jxyh_Jyxh')
        obs_temp = 0d0
        vec_l = -hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        vec_r = -hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        obs_temp = obs_temp + sum(vec_l * matmul(spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

!! ---------------------------spincurrent_hatree_only_correlation_end----------------------------!!

!! ---------------------------spin_momentum_correlation_begin----------------------------!!
        call spin_momentum_matrix(pc1,pc2,spinM_mat,4,hatree_only = .false.) ! spin momentum mat
      !! corf: den-den
        handle = this%get_handle('den_den')
        obs_temp = 0d0
        vec_l = [1,0,0,1]/2d0
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinM_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf: SC-SC
        handle = this%get_handle('SC_SC')
        obs_temp = 0d0
        obs_temp = obs_temp + (g_h(s2_u,s1_u) * g_h(s2_d,s1_d)- g_h(s2_d,s1_u)*g_h(s2_u,s1_d))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf6: Sz-Sz
        handle = this%get_handle('Sz_Sz')
        obs_temp = 0d0
        vec_l = [1,0,0,-1]/2d0
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinM_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf7: Sx-Sx
        handle = this%get_handle('Sx_Sx')
        obs_temp = 0d0
        vec_l = [0,1,1,0]/2d0
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinM_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf6: Sy-Sy
        handle = this%get_handle('Sy_Sy')
        obs_temp = 0d0
        vec_l = (0,1) * [0,-1,1,0]/2d0
        vec_r = vec_l
        obs_temp = obs_temp + sum(vec_l * matmul(spinM_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)
!! ---------------------------spin_momentum_correlation_end----------------------------!!

!! ---------------------------physical_spincurrent_correlation_begin-----------------------------!!
        if(allocated(vec_l)) deallocate(vec_l, vec_r)
        allocate(vec_l(8), vec_r(8))
      !! corf: physical y-bond spin-x current correlation
        call physical_spin_current_matrix(pc1,pc1_y,pc2,pc2_y,phy_spinJxy_mat,hatree_only = .false.) ! y-bond's spinJ mat
        handle = this%get_handle('cJyx_cJyx')
        obs_temp = 0d0
        vec_l(1:4) = hop * [(0,-1d0), (0,-1d0), (0,1d0), (0,1d0)] ! spin-x current vector
        vec_l(5:8) = ep_parameter * boson_field(bf1_y,time) * [1d0, 1d0, 1d0, 1d0] ! kinetic operator vector
        vec_r(1:4) = vec_l(1:4)
        vec_r(5:8) = ep_parameter * boson_field(bf2_y,time) * [1d0, 1d0, 1d0, 1d0] ! kinetic operator vector
        ! check whether pc1&pc1_y or pc2&pc2_y cross the boundary
        if(pc1%dpos(2)==Lat%dlength(2)-1) then
          vec_l(5:8) = vec_l(5:8) * [Lat%BC_phase(2), conjg(Lat%BC_phase(2)), Lat%BC_phase(2), conjg(Lat%BC_phase(2))]
        end if
        if(pc2%dpos(2)==Lat%dlength(2)-1) then
          vec_r(5:8) = vec_r(5:8) * [Lat%BC_phase(2), conjg(Lat%BC_phase(2)), Lat%BC_phase(2), conjg(Lat%BC_phase(2))]
        end if
        obs_temp = obs_temp + sum(vec_l * matmul(phy_spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      !! corf: physical x-bond spin-y current correlation
        call physical_spin_current_matrix(pc1,pc1_x,pc2,pc2_x,phy_spinJxy_mat,hatree_only = .false.) ! x-bond's spinJ mat
        handle = this%get_handle('cJxy_cJxy')
        obs_temp = 0d0
        vec_l(1:4) = hop * [-1d0, 1d0, 1d0, -1d0] ! spin-y current vector
        vec_l(5:8) = ep_parameter * boson_field(bf1_x,time) * [1d0, 1d0, 1d0, 1d0] ! kinetic operator vector
        vec_r(1:4) = vec_l(1:4)
        vec_r(5:8) = ep_parameter * boson_field(bf2_x,time) * [1d0, 1d0, 1d0, 1d0] ! kinetic operator vector
        ! check whether pc1&pc1_x or pc2&pc2_x cross the boundary
        if(pc1%dpos(1)==Lat%dlength(1)-1) then
          vec_l(5:8) = vec_l(5:8) * [Lat%BC_phase(1), conjg(Lat%BC_phase(1)), Lat%BC_phase(1), conjg(Lat%BC_phase(1))]
        end if
        if(pc2%dpos(1)==Lat%dlength(1)-1) then
          vec_r(5:8) = vec_r(5:8) * [Lat%BC_phase(1), conjg(Lat%BC_phase(1)), Lat%BC_phase(1), conjg(Lat%BC_phase(1))]
        end if
        obs_temp = obs_temp + sum(vec_l * matmul(phy_spinJxy_mat, vec_r))/N_cell
        call this%record_field_entry(handle, ind, obs_temp)

      end do !for c2
    end do !for c1

      !! END MEASUREMENT

  END SUBROUTINE ms_take_measurement

  subroutine ms_k_space_obs(this)
    implicit none
    class(MeasurementSystem), intent(inout) :: this
    integer :: i
    complex(dp) :: k_data(lat%N_cell)
    k_data = 0.0d0
    do i = 1, this%nobs
      if (this%kinds(i) == KIND_FIELD) then
        ! create k-space data
        k_data = matmul(lat%k_phase, this%data(this%off_lo(i):this%off_hi(i), this%cur_meas, this%cur_bin))/lat%N_cell
        ! record k-space data
        call this%record_field(this%get_handle(trim(this%names(i))//'_k'), k_data)
      end if
    end do
  end subroutine ms_k_space_obs

!======================================================================!
!                             NAMELIST INIT                             !
!======================================================================!

  subroutine ms_init_meas(this)
    implicit none
    class(MeasurementSystem), intent(inout) :: this
    this%nbins = nbin_per_core
    this%nmeas = nmeas_per_bin
    this%cur_bin = 1
    this%cur_meas = 1
    call this%init_from_namelist(nml_file)
  end subroutine ms_init_meas

  subroutine ms_init_from_namelist(this, filename)
    class(MeasurementSystem), intent(inout) :: this
    character(*),            intent(in)    :: filename
    integer :: lun, ios
    character(len=NAME_LEN) :: name, obskind
    integer :: cap
    namelist /observable/ name, obskind
    print*,this%nbins,this%nbins
    if (this%nbins<=0 .or. this%nmeas<=0) error stop "init_from_namelist: nbins/nmeas must be >0"
    this%field_width = product(Lat%dlength)
    ! 初始容量（动态增长，避免频繁重分配）
    cap = 50
    allocate(this%names(cap), this%kinds(cap), this%widths(cap), this%off_lo(cap), this%off_hi(cap))
    this%nobs = 0

    open(newunit=lun, file=filename, status='old', action='read',position = 'rewind' ,iostat=ios)
    if (ios/=0) error stop "init_from_namelist: cannot open file: "//trim(filename)

    do
      name = ''; obskind = ''
      read(lun, nml=observable, iostat=ios)
      if (ios /= 0) exit  ! 到末尾或无更多组
      call tolower_inplace(obskind)
      call trim_inplace(name)
      if (len_trim(name)==0) error stop "observable group missing name"
      select case (obskind)
      case ('scalar','s','sc','1')
        call push_obs(this, name, KIND_SCALAR)
      case ('field','f','array','a')
        if (this%field_width<=0) error stop "observable obskind=field but field_width not provided"
        call push_obs(this, name, KIND_FIELD)
        call push_obs(this, trim(name)//'_k', KIND_FIELD)
      case default
        error stop "unknown observable obskind: "//trim(obskind)//" (name="//trim(name)//")"
      end select
    end do
    close(lun)

    if (this%nobs==0) error stop "no observable groups found in namelist"
    print*,"MeasurementSystem: found ", this%nobs, " observables."
    call this%finalize_layout()
    call this%init_storage()
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
    print*,"MeasurementSystem: total data width = ", this%total_width
    print*,"MeasurementSystem: data layout:"
    do i=1,this%nobs
      !print*,"  ", trim(this%names(i)), " (kind=", this%kinds(i), ", width=", this%widths(i), &
      !& ", offset=[", this%off_lo(i), ":", this%off_hi(i), "])"
    end do
  end subroutine ms_finalize_layout

  subroutine ms_init_storage(this)
    class(MeasurementSystem), intent(inout) :: this
    if (.not.this%finalized) error stop "init_storage: call finalize_layout first"
    if (this%nbins<=0 .or. this%nmeas<=0) error stop "init_storage: nbins/nmeas must be >0"
    allocate(this%data(this%total_width, this%nmeas,this%nbins))
    this%data = 0.0_dp
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
    if (this%cur_bin > this%nbins) error stop "begin_measure: no more bins left"
    call this%take_measurement(time)
    ! next loop
    if ((forward .and. time == ntime/2 + meas_number) .or. ((.not. forward) .and. time == ntime/2 - meas_number)) then
      ! average over all sites and all measurements in one loop
      this%data(:,this%cur_meas, this%cur_bin) = this%data(:,this%cur_meas, this%cur_bin)/(2*meas_number + 1)
      call this%k_space_obs()
      !this%data(:,this%cur_meas, this%cur_bin) = this%data(:,this%cur_meas, this%cur_bin)/real(Lat%N_cell,dp)
      this%cur_meas = this%cur_meas + 1
      if (this%cur_meas > this%nmeas) then
        ! next bin
        this%cur_meas = 1
        this%cur_bin = this%cur_bin + 1
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
    h = -1
    return

  end function ms_get_handle

  subroutine ms_get_range_by_handle(this, handle, lo, hi, width, obskind)
    class(MeasurementSystem), intent(in)  :: this
    integer,                 intent(in)  :: handle
    integer,                 intent(out) :: lo, hi, width, obskind
    if (handle<1 .or. handle>this%nobs) error stop "get_range_by_handle: handle out of range"
    lo = this%off_lo(handle); hi = this%off_hi(handle)
    width = this%widths(handle); obskind  = this%kinds(handle)
  end subroutine ms_get_range_by_handle

  subroutine ms_record_scalar(this, handle, value, bin, meas)
    class(MeasurementSystem), intent(inout) :: this
    integer,                 intent(in)    :: handle
    complex(dp),                intent(in)    :: value
    integer,      intent(in), optional    :: bin, meas
    integer :: lo, hi, width, obskind, b, m
    if(handle == -1) return
    call this%get_range_by_handle(handle, lo, hi, width, obskind)
    if (obskind /= KIND_SCALAR) error stop "record_scalar: not a scalar handle"
    b = this%cur_bin
    m = this%cur_meas
    if (b<1 .or. b>this%nbins) error stop "record_scalar: bin out of range"
    if (m<1 .or. m>this%nmeas) error stop "record_scalar: meas out of range"
    this%data(lo, m, b) = this%data(lo, m, b) + value

  end subroutine ms_record_scalar

  subroutine ms_record_field(this, handle, arr, bin, meas)
    class(MeasurementSystem), intent(inout) :: this
    integer,                 intent(in)    :: handle
    complex(dp),                intent(in)    :: arr(:)   ! flattened field (length = field_width)
    integer,      intent(in), optional    :: bin, meas
    integer :: lo, hi, width, obskind, b, m
    if(handle == -1) return
    call this%get_range_by_handle(handle, lo, hi, width, obskind)
    if (obskind /= KIND_FIELD) error stop "record_field: not a field handle"
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
    complex(dp),                intent(in)    :: value
    integer,      intent(in), optional    :: bin, meas
    integer :: lo, hi, width, obskind, b, m
    if(handle == -1) return
    call this%get_range_by_handle(handle, lo, hi, width, obskind)
    if (obskind /= KIND_FIELD) error stop "record_field_entry: not a field handle"
    if (idx<1 .or. idx>width) error stop "record_field_entry: idx out of range"
    b = this%cur_bin
    m = this%cur_meas
    if (b<1 .or. b>this%nbins) error stop "record_field_entry: bin out of range"
    if (m<1 .or. m>this%nmeas) error stop "record_field_entry: meas out of range"
    this%data(lo+idx-1, m, b) = this%data(lo+idx-1, m, b) + value
  end subroutine ms_record_field_entry

!======================================================================!
!                               ANALYSIS                                !
!======================================================================!

  subroutine ms_data_analyse(this)
    class(MeasurementSystem), intent(inout) :: this
    ! 这里可以添加更多分析功能
    if (.not.allocated(this%data)) error stop "data_analyse: storage not initialized"
    if (.not.allocated(this%mean)) allocate(this%mean(this%total_width, this%nbins))
    call this%compute_all_bin_means()
    call this%output_array_data()
    call this%output_reading_guide()
  end subroutine ms_data_analyse

  subroutine ms_compute_all_bin_means(this)
    ! 计算所有观测量在每个 bin 的均值（沿第 3 维求平均）
    ! 返回: means(total_width, nbins)
    class(MeasurementSystem), intent(inout)  :: this
    integer :: i_b, i_m

    this%mean = 0.0_dp
    do i_b = 1, this%nbins
      do i_m = 1, this%nmeas
        this%mean(:, i_b) = this%mean(:, i_b) + real(this%data(:, i_m, i_b))
      end do
      this%mean(:, i_b) = this%mean(:, i_b) / real(this%nmeas, dp)
    end do
  end subroutine ms_compute_all_bin_means

  subroutine ms_output_array_data(this)

    class(MeasurementSystem), intent(in) :: this
    integer :: i_b,ios,line
    character(len = 100) :: ci1,ci2,str

    if (.not.allocated(this%mean)) error stop "output_array_data: mean not computed"
#ifdef MPI
    write (ci1, '(1i4)') myid
    str = 'out_core'//trim(adjustl(ci1))//'.dat'
#else
    str = 'out_core.dat'
#endif
    str = trim(adjustl(output_addr))//'data/'//trim(adjustl(str))
    write(ci2, '(1i4)') this%nbins
    open(unit=10, file=str, status='replace', iostat=ios)
    if ( ios /= 0) error stop "output_array_data: cannot open file: "//trim(str)
    do line = 1, this%total_width
      ! output array data, one index per line,do not specify observable name.
      ! The output data is (total_width, nbins)
      write (unit=10,FMT ='('//trim(adjustl(ci2))//output_format//')',iostat = ios) this%mean(line,:)
    end do
    close(10)
  end subroutine ms_output_array_data

  subroutine ms_output_reading_guide(this)
    class(MeasurementSystem), intent(in) :: this
    integer :: i, lun, ios
    character(len=100) :: ci1,ci2,str
#ifdef MPI
    if(myid/=0) return
#endif
    str = trim(adjustl(output_addr))//'data/'//reading_guide
    open(newunit=lun, file=str, status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop "output_reading_guide: cannot open file: "//trim(str)
    write(lun, *) '&basics'
    write(lun, *) 'n_obs=', this%nobs
    write(lun, *) 'total_width=', this%total_width
    write(lun, *) 'n_bins=', this%nbins
    write(lun, *) 'n_meas_per_bin=', this%nmeas
    write(lun, *) 'field_width=', this%field_width
#ifdef MPI
    write(lun, *) 'n_cores=', numprocs
    write(lun, *) 'MPI_nblock=',MPI_nblock
    write(lun, *) 'MPI_one_block=',MPI_one_block
#endif
    write(lun, *) 'output_format="', trim(adjustl(output_format))//'"'
    write(lun, *) '/'
    write(lun, *) ''

    do i=1,this%nobs
      write(lun,*) '&obs'
      write(lun,*) 'name="',trim(this%names(i))//'"'
      write(lun,*) 'kind=', this%kinds(i)
      write(lun,*) 'width=', this%widths(i)
      write(lun,*) 'offset_lo=', this%off_lo(i)
      write(lun,*) 'offset_hi=', this%off_hi(i)
      write(lun,*) '/'
      write(lun,*) ''
    end do
    close(lun)

  end subroutine ms_output_reading_guide
!======================================================================!
!                        USER MEASUREMENT PLACEHOLDER                   !
!======================================================================!
!
!======================================================================!

end module dqmc_measurements
