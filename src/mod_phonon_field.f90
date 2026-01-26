Module pf_setting
  !use input
  use nrtype
  use lattice
  use matrixlib

  ABSTRACT INTERFACE
    SUBROUTINE get_expKV_f(expKV, ppf,n_pla,bf_value,inv)
      IMPORT :: pf_type, dp
      TYPE(pf_type), INTENT(IN) :: ppf
      INTEGER, INTENT(IN) :: n_pla
      LOGICAL, INTENT(IN) :: inv
      REAL(dp),INTENT(IN) :: bf_value
      COMPLEX(dp), INTENT(OUT) :: expKV(ppf%dim, ppf%dim)
    END SUBROUTINE get_expKV_f
  END INTERFACE

  ! 2. 声明过程指针
  PROCEDURE(get_expKV_f), POINTER :: get_expKV => NULL()
contains

  subroutine set_pf_main()
    implicit none
    integer :: i_pf,i_pla,i_site,i_site2
    type(pf_type),pointer :: ppf

    if(.not.allocated(pf_list)) allocate(pf_list(n_phonon_field),ppf_list(n_phonon_field),pf_data_list(n_phonon_field))
    call setup_expKV(model_name)
    call readin_pf_basic()
    !print*,'Finished reading in pf basic info.'
    call set_pf_detail()

    ! debug
    if(.false.) then
      ! print the pf_site_list information
      do i_pf = 1, 1
        ppf => pf_list(i_pf)
        print*,'phonon field id:',ppf%id,'dim:',ppf%dim,'n_plaquette:',ppf%n_plaquette
        do i_pla = 1, ppf%n_plaquette
          write(*,'(A,I4,A)',advance='no') 'plaquette ',i_pla,' sites:'
          do i_site = 1, ppf%dim
            write(*,'(I6)',advance='no') ppf%p_data%pla_site_list(i_site,i_pla)
          end do
          print*,''
        end do
        write(*,'(A,1f18.6)',advance='no')'bf:',boson_field(ppf%p_data%bf_list(1),2)*ppf%V_coe
        do i_site = 1,ppf%dim
          print*,' '
          do i_site2 = 1,ppf%dim
            write(*,'(A,2f10.6,A)',advance='no')'  (',ppf%p_data%expKV(i_site,i_site2,1,2),')  '
          end do

        end do
      end do

    end if
  end subroutine set_pf_main

  subroutine readin_pf_basic()
    implicit none
    integer :: i_pf
    integer :: n_plaquette
    integer :: nml_unit  ! 定义一个变量来存设备号
    integer,dimension(n_phonon_field) :: id, pf_dim
    real(dp),dimension(n_phonon_field) ::  n_cover
    logical,dimension(n_phonon_field) :: K_exist, V_exist
    integer :: ios
    type(pf_type),pointer :: ppf
    namelist /pf_basic/ id, K_exist, V_exist, pf_dim, n_cover

    id = 0
    K_exist = .false.
    V_exist = .false.
    pf_dim = 0
    n_cover = 0
    open(newunit=nml_unit, file=trim(nml_file), status='old', &
         action='read', iostat=ios)
    if (ios /= 0) then
      print *, "Error opening namelist file:", trim(nml_file)
      print*,'ios:',ios
      stop
    end if
    read(nml_unit, nml=pf_basic,iostat=ios)
    if (ios /= 0) then
      print *, "Error reading pf_basic namelist"
      print*,'ios:',ios
      stop
    end if
    !print*,'Readin phonon field basic info:', id, K_exist, V_exist, pf_dim, n_cover

    do i_pf = 1, n_phonon_field

      ppf => pf_list(i_pf)
      !rewind(10)
      n_plaquette = int(Lat%Ns * n_cover(i_pf)/pf_dim(i_pf))
      ! 把读到的临时变量赋值到 pf
      !print*,'i_pf:',i_pf,'id:',id(i_pf),'K_exist',K_exist(i_pf),'V_exist',V_exist(i_pf), &
      !  'dim:',pf_dim(i_pf),'n_plaquette:',n_plaquette
      ppf%id = id(i_pf)
      ppf%K_exist = K_exist(i_pf)
      ppf%V_exist = V_exist(i_pf)
      ppf%dim = pf_dim(i_pf)
      ppf%n_plaquette = n_plaquette
      !ppf%site_list(1:n_plaquette, 1:n_phonon_field) = reshape(site_list, &
      !    shape=[n_plaquette, n_phonon_field])

    end do

    close(nml_unit)
    ppf_list => pf_list

  end subroutine readin_pf_basic

  subroutine set_pf_detail()
    implicit none
    integer :: i_pf,bf_index
    integer :: nml_unit
    character(len=16) :: group_name
    complex(dp),allocatable :: Vmatrix(:),Kmatrix(:)
    integer, allocatable ::  pla_tsl_dvec_uc(:), pla_offset_dvec_uc(:)
    integer,allocatable :: pla_int_subsites(:)
    integer :: ios
    type(pf_type),pointer :: ppf
    complex(dp) :: K_coe,V_coe
    namelist /pf_pos/&
      pla_tsl_dvec_uc, pla_offset_dvec_uc, pla_int_subsites
    namelist /pf_K/&
      K_coe,Kmatrix
    namelist /pf_V/&
      bf_index,V_coe,Vmatrix
    open(newunit=nml_unit, file=trim(nml_file), status='old',action='read',iostat=ios)
    if (ios /= 0) then
      print *, "Error opening namelist file:", trim(nml_file)
      print*,'ios:',ios
      stop
    end if
    !print*,'Setting phonon field detail info...'
    do i_pf = 1, n_phonon_field

      write(group_name, '("pfdetail", I0)') i_pf
      ppf => pf_list(i_pf)
      allocate(Vmatrix(ppf%dim*ppf%dim), &
               Kmatrix(ppf%dim*ppf%dim), &
               pla_tsl_dvec_uc(Lat%dim*Lat%dim), &
               pla_offset_dvec_uc(Lat%dim), &
               pla_int_subsites((Lat%dim+1)*ppf%dim))
      Vmatrix = 0.0d0
      Kmatrix = 0.0d0
      pla_tsl_dvec_uc = 0
      pla_offset_dvec_uc = 0
      pla_int_subsites = 0
      bf_index = 0
      K_coe = 0.0d0
      V_coe = 0.0d0
      read(nml_unit, nml=pf_pos, iostat=ios)
      if (ios /= 0) then
        print *, "Error reading pf_pos namelist for phonon field ", i_pf
        print*,'ios:',ios
        stop
      end if
      if(ppf%K_exist) then
        read(nml_unit, nml=pf_K, iostat=ios)
        if (ios /= 0) then
          print *, "Error reading pf_K namelist for phonon field ", i_pf
          print*,'ios:',ios
          stop
        end if
        ! K_coe = K_coe * hop
      end if
      if(ppf%V_exist) then
        read(nml_unit, nml=pf_V, iostat=ios)
        if (ios /= 0) then
          print *, "Error reading pf_V namelist for phonon field ", i_pf
          print*,'ios:',ios
          stop
        end if
        ! V_coe = ep_parameter * V_coe
      end if
      ppf%p_data => pf_data_list(ppf%id)
      call set_typed_pf(ppf,ppf%dim,K_coe,V_coe, &
                        Kmatrix,Vmatrix,pla_tsl_dvec_uc, pla_offset_dvec_uc, pla_int_subsites, bf_index)
      !allocate(ppf%Vmatrix,mold = Vmatrix)
      ! 把读到的临时变量赋值到 pf

      deallocate(Vmatrix, Kmatrix, &
                 pla_tsl_dvec_uc, pla_offset_dvec_uc, pla_int_subsites)

      !ppf%site_list(1:n_plaquette, 1:n_phonon_field) = reshape(site_list, &
      !    shape=[n_plaquette, n_phonon_field])

    end do

    close(nml_unit)
    ppf_list => pf_list

  end subroutine set_pf_detail

  subroutine readin_slater_pf()
    implicit none
    integer ::  pf_dim ,ios
    integer :: nml_unit
    real(dp) :: n_cover
    complex(dp),allocatable :: Vmatrix(:),Kmatrix(:)
    integer, allocatable ::  pla_tsl_dvec_uc(:), pla_offset_dvec_uc(:), pla_int_subsites(:)
    complex(dp) :: K_coe,V_coe
    type(pf_type),pointer :: ppf
    namelist /slater/ pf_dim, n_cover
    namelist /slater_detail/ K_coe, Kmatrix, pla_tsl_dvec_uc, pla_offset_dvec_uc, pla_int_subsites
    ! this subroutine is used to set the Slater states for trial wavefunciton
    ! slater field is read in from the namelist_slater.nml
    ! usually it is one matrix(plaquette) of size (Ns,Ns)
    ! first readin slater and slater_detail, the data is in namelist_pf.nml
    open(newunit=nml_unit, file=trim(nml_file), status='old',action='read',iostat=ios)
    read(nml_unit, nml=slater, iostat=ios)
    if (ios /= 0) then
      print *, "Error reading slater namelist"
      print*,'ios:',ios
      stop
    end if
    ppf => slater_pf
    ppf%dim = pf_dim
    ppf%n_plaquette = int(Lat%Ns * n_cover / pf_dim)
    ppf%id = -1 ! set the id to 1 for slater
    ppf%K_exist = .true.
    ppf%V_exist = .false.
    allocate(Vmatrix(pf_dim*pf_dim), &
             Kmatrix(pf_dim*pf_dim), &
             pla_tsl_dvec_uc(Lat%dim*Lat%dim), &
             pla_offset_dvec_uc(Lat%dim), &
             pla_int_subsites((Lat%dim+1)*pf_dim))
    K_coe = 0.0d0
    V_coe = 0.0d0
    Vmatrix = 0.0d0
    Kmatrix = 0.0d0
    pla_tsl_dvec_uc = 0
    pla_offset_dvec_uc = 0
    pla_int_subsites = 0
    read(nml_unit, nml=slater_detail, iostat=ios)
    if (ios /= 0) then
      print *, "Error reading slater_detail namelist"
      print*,'ios:',ios
      stop
    end if
    ppf%p_data => slater_pf_data
    call set_typed_pf(ppf, pf_dim,K_coe,V_coe,Kmatrix, Vmatrix, pla_tsl_dvec_uc, pla_offset_dvec_uc, pla_int_subsites)
    close(nml_unit)
  end subroutine readin_slater_pf

  subroutine set_typed_pf(ppf,dim,K_coe,V_coe,Kmatrix,Vmatrix, pla_tsl_dvec_uc, pla_offset_dvec_uc, pla_int_subsites,bf_index)
    implicit none
    ! this subroutine sets ppf%Vmatrix, Kmatrix, pla_tsl_dvec_uc, pla_offset_dvec_uc, pla_int_subsites
    ! and also the corresponding pf_data
    ! the ppf's dim and n_plaquette should be set before calling this subroutine and dim = ppf%dim
    ! the input Vmatrix, Kmatrix, pla_tsl_dvec_uc, pla_offset_dvec_uc, pla_int_subsites should be in 1d array
    ! the reshape function is used to convert them to 2d matrix
    ! pla_tsl_dvec_uc contains the translational dvec of plaquettes in unit cell basis, shape (Lat%dim, Lat%dim)
    ! pla_offset_dvec_uc is in shape (Lat%dim), contains the offset dvec of plaquettes in unit cell basis
    ! the above two only gives the dpos of the representative site of each plaquette, also defined as the first site in plaquette
    ! pla_int_subsites is in shape (Lat%dim+1, ppf%dim), contains ppf%dim of lat%dim+1 vectors
    ! each vector contains (the relative dpos between UCs of each subsite and the first site + the subsite index in uc)
    ! the first site's vector is always (0,0,0,..., subsite_index)
    ! bf_index is optional, if present, the bf_list and expKV will be set
    ! cross_boundary and BC_phase are also set in this subroutine
    ! for each plaquette, we check whether any directional pair of sites cross the boundary, if yes, set cross_boundary to true
    ! and the BC_phase matrix is set accordingly
    type(pf_type),intent(inout) :: ppf
    integer,intent(in) :: dim
    type(pf_data_type),pointer :: p_data_temp
    integer,intent(in),optional :: bf_index
    complex(dp),intent(in) :: K_coe,V_coe,Vmatrix(dim*dim),Kmatrix(dim*dim)
    integer, intent(in) ::pla_tsl_dvec_uc(Lat%dim*Lat%dim),pla_offset_dvec_uc(Lat%dim)
    integer, intent(in) :: pla_int_subsites((Lat%dim+1)*ppf%dim)
    integer :: i_site,j_site,i_cell, i_plaquette,pla,pla_count,cell_index,i,j,temp_id,k
    integer :: pla_pos_vec(lat%dim),cell_dpos(lat%dim),i_dpos(lat%dim),j_dpos(lat%dim)
    real(dp) :: inv_tsl_vec(lat%dim,lat%dim),cell_xpos(lat%dim)
    allocate(ppf%Vmatrix(dim,dim), &
             ppf%Kmatrix(dim,dim))
    allocate(ppf%pla_tsl_dvec_uc(Lat%dim,Lat%dim), &
             ppf%pla_offset_dvec_uc(Lat%dim), &
             ppf%pla_int_subsites(Lat%dim+1,dim))
    ppf%K_coe = K_coe
    ppf%V_coe = V_coe
    ppf%Vmatrix = reshape(Vmatrix, shape=[dim, dim])
    ppf%Kmatrix = reshape(Kmatrix, shape=[dim, dim])
    ppf%pla_tsl_dvec_uc = reshape(pla_tsl_dvec_uc, shape=[Lat%dim, Lat%dim])
    ppf%pla_offset_dvec_uc = pla_offset_dvec_uc
    ppf%pla_int_subsites = reshape(pla_int_subsites, shape=[Lat%dim+1, dim])
    if (.not. associated(ppf%p_data)) then
      allocate(ppf%p_data)
    end if
    p_data_temp => ppf%p_data
    ! allocate the p_data
    allocate(p_data_temp%pla_site_list(ppf%dim,ppf%n_plaquette),p_data_temp%boundary_crossing(ppf%n_plaquette), &
             p_data_temp%BC_phases(ppf%dim,ppf%dim,ppf%n_plaquette))
    p_data_temp%pla_site_list = -1
    p_data_temp%boundary_crossing = .false.
    p_data_temp%BC_phases = (1.0d0,0.0d0)
    ! run over all cells in the lattice and check if the cell is in one plaquette of the phonon field
    ! to check if the cell is in the plaquette, we need to check if the cell's dpos(denoted by v) is integer multiple of
    ! the pf translational dvec(denoted by M). That's to say M x = v has integer solution x.
    ! x = inv(M) * v, and see if mod(x,1) == 0.

    inv_tsl_vec = real(ppf%pla_tsl_dvec_uc)
    call inverse(Lat%dim, inv_tsl_vec)

    pla_count = 0
    do i_cell = 1, Lat%N_cell
      cell_dpos = p_cells(i_cell)%dpos - ppf%pla_offset_dvec_uc
      ! correction if the site cross the boundary
      do j = 1, Lat%dim
        if(cell_dpos(j) < -0.0001d0) then
          cell_dpos(j) = cell_dpos(j) + Lat%dlength(j)
        end if
        if(cell_dpos(j) >= Lat%dlength(j)-0.0001d0) then
          cell_dpos(j) = cell_dpos(j) - Lat%dlength(j)
        end if
      end do
      ! calculate the position of the site in the translational vector
      cell_xpos = matmul(inv_tsl_vec,cell_dpos)
      if(any(abs(nint(cell_xpos) - cell_xpos) > 1e-6)) cycle ! if the site is not in the pf, skip it
      pla_count = pla_count + 1
      ! --- 新增防越界检查 ---
      if (pla_count > ppf%n_plaquette) then
        print *, "FATAL ERROR in set_typed_pf: pla_count Exceeds limit!"
        print *, "pla_count =", pla_count, " limit (n_plaquette) =", ppf%n_plaquette
        stop
      end if
      ! --------------------
      ! record the subsite index of the first site in the plaquette(the representative site)
      p_data_temp%pla_site_list(1,pla_count) = p_cells(i_cell)%sites(ppf%pla_int_subsites(lat%dim+1,1))
    end do

    if(ppf%dim > 1) then
      ! loop over each plaquette and determine the plaquette index and pos vec first
      do i_plaquette = 1, ppf%n_plaquette
        pla_pos_vec = p_sites(p_data_temp%pla_site_list(1,i_plaquette))%uc_dpos
        do i_site = 2, ppf%dim
          ! calculate the position of unit cell that occupies the i-th site in the plaquette
          cell_dpos = pla_pos_vec + ppf%pla_int_subsites(1:lat%dim,i_site)
          ! check if the bond from site 1 to i_site cross the boundary
          do j = 1, Lat%dim

            if(cell_dpos(j) > Lat%dlength(j)-0.0001d0) then ! site cross positive boundary
              cell_dpos(j) = cell_dpos(j) - Lat%dlength(j)
            end if
            if(cell_dpos(j) < -0.0001d0) then ! site cross negative boundary
              cell_dpos(j) = cell_dpos(j) + Lat%dlength(j)
            end if
            !cell_dpos(j) = modulo(cell_dpos(j), Lat%dlength(j))
          end do
          ! calculate the index of the site
          call get_uc_index_from_dpos(cell_dpos,cell_index)
          ! store the site index in the pla_site_list
          p_data_temp%pla_site_list(i_site,i_plaquette) = p_cells(cell_index)%sites(ppf%pla_int_subsites(lat%dim+1,i_site))
        end do
      end do

      ! apply BC_phase if crossing boundary
      do i_plaquette = 1, ppf%n_plaquette
      do i = 1,ppf%dim
        do j = 1,ppf%dim
          i_site = p_data_temp%pla_site_list(i,i_plaquette)
          i_dpos = p_sites(i_site)%uc_dpos
          j_dpos = i_dpos - ppf%pla_int_subsites(1:lat%dim,i) + ppf%pla_int_subsites(1:lat%dim,j) ! relative dpos inferred
          ! check boundary crossing
          ! different subsite, no BC phase
          ! if(abs(ppf%pla_int_subsites(lat%dim+1,i) - ppf%pla_int_subsites(lat%dim+1,j)) < 0.01d0) cycle
          do k  = 1, Lat%dim
            if(j_dpos(k) >= Lat%dlength(k)-0.0001d0) then
              p_data_temp%boundary_crossing(i_plaquette) = .true.
              p_data_temp%BC_phases(i,j,i_plaquette) = p_data_temp%BC_phases(i,j,i_plaquette) * lat%BC_phase(k)
            end if
            if(j_dpos(k) < -0.0001d0) then
              p_data_temp%boundary_crossing(i_plaquette) = .true.
              p_data_temp%BC_phases(i,j,i_plaquette) = p_data_temp%BC_phases(i,j,i_plaquette) * conjg(lat%BC_phase(k))
            end if
          end do

        end do
      end do
      end do

      ! debug
      if(.false.) then
        do i_plaquette = 1, ppf%n_plaquette
          write(*,'(A,I4,A)',advance='no') 'plaquette ',i_plaquette,' sites:'
          do i_site = 1, ppf%dim
            write(*,'(I6)',advance='no') p_data_temp%pla_site_list(i_site,i_plaquette)
          end do
          print*,''
          print*,'boundary crossing:',p_data_temp%boundary_crossing(i_plaquette)
          print*,'BC phases matrix:'
          do i_site = 1, ppf%dim
            write(*,'(A)',advance='no') '  '
            do j = 1, ppf%dim
              write(*,'(A,2f10.6,A)',advance='no')'  (',p_data_temp%BC_phases(i_site,j,i_plaquette),')  '
            end do
            print*,''
          end do

        end do
        stop
      end if

    end if

    ! set the bf_list and expKV if bf_index is present
    if(present(bf_index)) then
      ! set the bf_list

      allocate(p_data_temp%bf_list(ppf%n_plaquette))
      do i_plaquette = 1, ppf%n_plaquette
        temp_id = p_sites(ppf%p_data%pla_site_list(1,i_plaquette))%uc_id
        p_data_temp%bf_list(i_plaquette) = (bf_index) + (temp_id - 1) * bf_sets
      end do
      ! set the expKV
      allocate(p_data_temp%expKV(ppf%dim,ppf%dim,ppf%n_plaquette,ntime),&
      & p_data_temp%expKV_inv(ppf%dim,ppf%dim,ppf%n_plaquette,ntime))
      call init_expKV(ppf)

    end if
  end subroutine set_typed_pf

  subroutine init_expKV(pf)
    implicit none
    type(pf_type),intent(inout) :: pf
    type(pf_data_type),pointer :: p_data_temp
    integer :: i_plaquette,p
    ! --- 调试代码开始 ---

    !if (.not. associated(pf%p_data)) then
    !  print *, "CRITICAL ERROR: ppf%p_data is NOT associated in init_expKV!"
    !  stop
    !end if
    !print *, "DEBUG: ppf%p_data is associated."
    ! --- 调试代码结束 ---
    p_data_temp => pf%p_data
    do p = 1, ntime
      do i_plaquette = 1, pf%n_plaquette
        if(pf%V_exist) then
          call get_expKV(p_data_temp%expKV(:,:,i_plaquette,p),pf,i_plaquette,&
          & boson_field(p_data_temp%bf_list(i_plaquette),p), .false.)
          call get_expKV(p_data_temp%expKV_inv(:,:,i_plaquette,p),pf,i_plaquette,&
          & boson_field(p_data_temp%bf_list(i_plaquette),p),.true.)
        else
          call get_expKV(p_data_temp%expKV(:,:,i_plaquette,p),pf,i_plaquette,0d0,.false.)
          call get_expKV(p_data_temp%expKV_inv(:,:,i_plaquette,p),pf,i_plaquette,0d0,.TRUE.)
        end if
      end do
    end do
    BLOCK ! debug block
      complex(dp) :: expKV1(pf%dim, pf%dim), expKV2(pf%dim, pf%dim)
      if(.false.) then
        ! compare the expKV_general and expKV_epsoc for epsoc model
        ! they should be the same
        do p = 1, ntime
        do i_plaquette = 1, pf%n_plaquette
          call get_expKV_general(expKV1, pf,i_plaquette, &
          & boson_field(p_data_temp%bf_list(i_plaquette),1), .false.)
          call get_expKV_epsoc(expKV2, pf,i_plaquette, &
          & boson_field(p_data_temp%bf_list(i_plaquette),1), .false.)
          if(maxval(abs(expKV1 - expKV2)) > 1e-6) then
            print*,'expKV mismatch for plaquette ',i_plaquette,'at time',p
            print*,'max difference:',maxval(abs(expKV1 - expKV2))
          end if
        end do
        END DO
      END IF
    END BLOCK
  end subroutine init_expkv

  SUBROUTINE setup_expKV(model_name)
    CHARACTER(LEN=*), INTENT(IN) :: model_name

    SELECT CASE (TRIM(model_name))
    CASE ('EPSOCZ')
      ! 指针指向优化版
      get_expKV => get_expKV_epsoc
      if(myid == 0) then
        print*,'Using optimized expKV for EPSOCZ model.'
      end if
    CASE DEFAULT
      ! 指针指向标准版
      get_expKV => get_expKV_general
      if(myid == 0) then
        print*,'Using general expKV calculation.'
      end if

    END SELECT
  END SUBROUTINE setup_expKV
  subroutine get_expKV_general(expKV, ppf,n_pla,bf_value,inv)
    implicit none
    type(pf_type), intent(in) :: ppf
    integer, intent(in) :: n_pla
    logical, intent(in) :: inv
    real(dp),intent(in) :: bf_value
    complex(dp), intent(out) :: expKV(ppf%dim, ppf%dim)

    !p_data => ppf%p_data
    !bf_index = p_data%bf_list(i_pla)
    expKV = 0d0
    if(ppf%V_exist .and. ppf%K_exist) then
      ! expKV = exp(-delt * (-t * K + gX * V))
      expKV = ppf%K_coe * ppf%Kmatrix + ppf%V_coe * bf_value * ppf%Vmatrix
    else if(ppf%V_exist .and. .not.ppf%K_exist) then
      ! expKV = exp(-delt * (gX * V))
      expKV = ppf%V_coe * bf_value * ppf%Vmatrix
    else if(.not.ppf%V_exist .and. ppf%K_exist) then
      ! expKV = exp(-delt * (-t * K))
      expKV = ppf%K_coe * ppf%Kmatrix
    else
      expKV = 0d0
    end if

    if(ppf%p_data%boundary_crossing(n_pla)) then
      expKV = expKV * ppf%p_data%BC_phases(:,:,n_pla)
    end if

    if(.not.inv) then
      expKV = -delt * expKV
      call expm(expKV,ppf%dim)
    else
      expKV = delt * expKV
      call expm(expKV,ppf%dim)
    end if
  end subroutine get_expKV_general

  subroutine get_expKV_epsoc(expKV, ppf,n_pla,bf_value,inv)
    ! this subroutine is used to get the expKV for epsoc model, which should be faster than the above one
    implicit none
    type(pf_type), intent(in) :: ppf
    integer, intent(in) :: n_pla
    logical, intent(in) :: inv
    real(dp),intent(in) :: bf_value
    complex(dp), intent(out) :: expKV(ppf%dim, ppf%dim)
    complex(dp) :: sub_h(2,2)
    real(dp) :: norm
    integer :: i,j
    expKV = 0d0
    ! for epsoc model, K_exist is always true, V_exist is always true
    ! exp(-delt * H) = exp(-delt * (-t * K + gX * V)) = cosh(delt * sqrt(t^2 + (gX * bf_value)^2)) * I
    ! - H * sinh(delt * sqrt(t^2 + (gX * bf_value)^2)) / sqrt(t^2 + (gX * bf_value)^2)
    ! or to say: exp(-delt * H) = cosh(norm) * I - delt * H * sinh(norm) / norm
    ! where norm = delt * sqrt(t^2 + (gX * bf_value)^2)
    ! sub_h is 2x2 matrix where the off-diagonal elements are non-zero
    norm = delt * sqrt((abs(ppf%K_coe))**2 + (abs(ppf%V_coe) * bf_value)**2)
    sub_h =- delt*(ppf%K_coe * ppf%Kmatrix(1:2,3:4) + ppf%V_coe * bf_value * ppf%Vmatrix(1:2,3:4))
    if(ppf%p_data%boundary_crossing(n_pla)) then
      sub_h = sub_h * ppf%p_data%BC_phases(1:2,3:4,n_pla)
    end if
    do i = 1, ppf%dim
      ! set the diagonal elements
      expKV(i,i) = cosh(norm)
    end do

    if(.not.inv) then
      expKV(1:2,3:4) = sub_h * sinh(norm) / norm
    else
      expKV(1:2,3:4) = -sub_h * sinh(norm) / norm
    end if
    expKV(3:4,1:2) = conjg(transpose(expKV(1:2,3:4)))

  end subroutine get_expKV_epsoc
end module pf_setting
