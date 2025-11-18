Module lattice
 use input
 use matrixlib
contains
 
  subroutine set_lattice_main()
    implicit none

    call readin_lattice_para()
    call setup_lattice()
    ! debug

  end subroutine set_lattice_main

  subroutine readin_lattice_para()

    implicit none
    integer :: i_pf,i_dim
    integer :: dim,n_subsite_uc
    integer, allocatable :: lat_size(:)
    real(8), allocatable :: lat_vec(:), site_vec(:)
    integer :: ios
    namelist /lat_dim/ dim,n_subsite_uc
    namelist /lat_spatial/  lat_size,lat_vec
    namelist /lat_subsite_vector/ site_vec
    !! Read lattice dimension and number of sites in unit cell
    dim = 0
    open (10, file=nml_file, status='old', position='rewind')
    !rewind(10)
    read (10, nml=lat_dim, iostat=ios)
    if (ios /= 0) then
      print *, "Error reading lattice dimension"
      print*,"ios = ", ios
      stop
    end if
    !print *, "Lattice dimension:", dim
    !! Read lattice vectors
    allocate(lat_size(dim),lat_vec(dim * dim),site_vec(n_subsite_uc * dim))
    rewind(10)
    read (10, nml=lat_spatial, iostat=ios)
    if (ios /= 0) then
      print *, "Error reading lattice spatial size and vectors"
      print*,"ios = ", ios
      stop 
    end if
    rewind(10)
    read (10, nml=lat_subsite_vector, iostat=ios)
    if (ios /= 0) then
      print *, "Error reading lattice unit cell vectors"
      print*,"ios = ", ios
      stop
    end if
    close(10)
    ! Initialize lattice type 
    lat%dim = dim
    !print *, "Lattice dimension from file:", lat%dim
    lat%n_subsite_uc = n_subsite_uc
    allocate(lat%dlength(dim), lat%rlength(dim))         
    allocate(lat%tsl_rvec(dim,dim))
    allocate(lat%subsite_rvec(dim,n_subsite_uc))
    lat%dlength = lat_size
    ! lat%tsl_rvec is the lattice vectors in spatial vectors
    lat%tsl_rvec = reshape(lat_vec, [dim, dim])
    lat%subsite_rvec = reshape(site_vec, [dim, n_subsite_uc])
    lat%N_cell = product(lat%dlength)
    lat%Ns = lat%N_cell * lat%n_subsite_uc
    lat%periodic = .true. ! Assuming periodic boundary conditions for now
    ! lat%rlength is the lat%dlength * norm of the lattice vectors in each dimension
    lat%rlength = real(lat%dlength, kind=8) * sqrt(sum(lat%tsl_rvec**2, dim=1)) ! assuming lat_vec is in the same order as dlength
    ! Debug output
    !print *, "Lattice dimension:", lat%dim
    !print *, "Lattice lengths:", lat%dlength
    !print *, "Lattice spatial vectors:"
    do i_dim = 1, dim
      print *, lat%tsl_rvec(:,i_dim)  
    end do 
    !print *, "Number of sites in unit cell:", lat%n_subsite_uc
    !print *, "site unit cell vectors:"
    do i_dim = 1, lat%n_subsite_uc
      print *, lat%subsite_rvec(:,i_dim)
    end do
  end subroutine readin_lattice_para

  subroutine setup_lattice()
    implicit none
    integer :: i, j, site_id
    ! This subroutine can be used to set up the lattice based on the parameters read in
    allocate(typed_sites(lat%Ns))
    allocate(typed_cells(lat%N_cell))
    do i = 1, lat%N_cell
      typed_cells(i)%id = i
      allocate(typed_cells(i)%dpos(lat%dim), typed_cells(i)%rpos(lat%dim), typed_cells(i)%sites(lat%n_subsite_uc))
      typed_cells(i)%dpos = 0
      typed_cells(i)%rpos = 0.0
      ! Calculate the position in lattice vectors
      ! i = 1+ (x1-1) + L1*((x2-1) + L2*((x3-1) + ...)), now we need to find x1, x2, ... 
      call get_dpos_from_uc_index(i,lat%dim,lat%dlength, typed_cells(i)%dpos)
      typed_cells(i)%rpos = reshape(matmul(lat%tsl_rvec, &
          reshape(typed_cells(i)%dpos, [lat%dim, 1])), [lat%dim])
      ! Now we have the position of the cell, we can set the sites in this cell
        do j = 1, lat%n_subsite_uc
            site_id = (i - 1) * lat%n_subsite_uc + j
            typed_sites(site_id)%id = site_id
            typed_sites(site_id)%uc_id = i
            typed_sites(site_id)%subsite_id = j
            allocate(typed_sites(site_id)%uc_dpos(lat%dim), &
                typed_sites(site_id)%uc_rpos(lat%dim), &
                typed_sites(site_id)%site_rpos(lat%dim))
            typed_sites(site_id)%uc_dpos = typed_cells(i)%dpos
            typed_sites(site_id)%uc_rpos = typed_cells(i)%rpos
            typed_sites(site_id)%site_rpos = typed_cells(i)%rpos+ &
                lat%subsite_rvec(:,j) ! position vec relative to the cell position
            typed_sites(site_id)%p_uc => typed_cells(i) ! pointer to the cell

        end do 
      typed_cells(i)%sites = [ ( j, j=1, lat%n_subsite_uc ) ] + (i-1)*lat%n_subsite_uc

    end do

    ! calculate the reciprocal lattice vectors
    allocate(lat%recip_vec(lat%dim,lat%dim),lat%FBZ(lat%dim,lat%N_cell),lat%k_phase(lat%N_cell,lat%N_cell))
    call calculate_1stBZ(lat%tsl_rvec, lat%recip_vec, lat%dim,lat%FBZ, lat%k_phase)

    ! Assign the pointer to the typed_sites
    allocate(p_sites(size(typed_sites)))
    allocate(p_cells(size(typed_cells)))
    p_sites => typed_sites
    p_cells => typed_cells
    ! print lattice information for debug
    !print *, "Total number of sites in the lattice:", lat%Ns
    !print *, "Lattice periodicity:", lat%periodic
    !print *, "Lattice unit cell count:", lat%N_cell
    !print *, "Lattice size:", lat%dlength
    !print *, "Lattice spatial lengths:", lat%rlength
    !print *, "Lattice spatial vectors:"
    do i = 1, lat%dim
     ! print *, lat%tsl_rvec(:,i)
    end do
    !print *, "Lattice setup completed."
  end subroutine setup_lattice

  subroutine calculate_1stBZ(rvec,kvec,n,FBZ,k_phase)
    implicit none
    integer,intent(in) :: n
    real(8) :: rvec(n,n),kvec(n,n),FBZ(n,lat%N_cell)
    complex(8) :: k_phase(lat%N_cell,lat%N_cell)
    integer :: i_dim,i_k,i_cell,k_dpos(n)
    real(8) :: k_offset(n),del_kvec(n,n)
    ! rvec(:,i) is the i-th lattice vector in spatial vectors
    ! kvec(:,i) is the i-th lattice vector in k-space vectors
    ! MATMUL(kvec(:,i)^T,  rvec(:,j)) = 2*pi*delta_ij
    ! kvec = 2pi * inv(rvec^T)
    ! n is the dimension of the lattice
    kvec = transpose(rvec)
    call inverse(n,kvec)
    kvec = kvec * 2.0d0 * acos(-1.0d0)
    ! ensure all the kvecs point in the same general direction 
    if(dot_product(kvec(:,1),rvec(:,1))/(norm2(kvec(:,1))*norm2(rvec(:,1))) < -0.00001d0) kvec(:,1) = -kvec(:,1)
    do i_dim = 1,n
      if(dot_product(kvec(:,i_dim),kvec(:,1))/(norm2(kvec(:,i_dim))*norm2(kvec(:,1))) < -0.00001d0) then
        kvec(:,i_dim) = -kvec(:,i_dim)
      end if
    end do
    ! calculate the first Brillouin zone
    k_offset = 0d0 ! k_offset is introduced to deal with odd size lattice
    del_kvec = 0d0
    do i_dim = 1, n
      k_offset = k_offset + mod(lat%dlength(i_dim),2) * kvec(:,i_dim)/lat%dlength(i_dim)/2.0d0
      k_offset = k_offset - 0.5d0 * kvec(:,i_dim) ! shift the FBZ to be centered at Gamma point
      del_kvec(:,i_dim) = kvec(:,i_dim)/lat%dlength(i_dim)
    end do
    do i_k = 1, lat%N_cell
        call get_dpos_from_uc_index(i_k,n,lat%dlength,k_dpos)
        FBZ(:,i_k) = matmul(del_kvec,k_dpos) + k_offset
    end do
    ! calculate the k_phase matrix
    k_phase = 0d0
    do i_k = 1, lat%N_cell
      do i_cell = 1, lat%N_cell
        k_phase(i_k,i_cell) = exp(complex(0d0,1d0)*dot_product(FBZ(:,i_k),typed_cells(i_cell)%rpos))
      end do
    end do
  end subroutine

  function get_site_index_from_rpos(rpos) result(re)
    implicit none
    real(8), intent(in) :: rpos(:)
    integer :: i, site_index,re
    real(8) :: diff
    ! This function returns the site index from the position in spatial vectors
    site_index = -1
    do i = 1, size(p_sites)
      diff = norm2(p_sites(i)%site_rpos - rpos)
      if (diff < 1.0e-6) then
        site_index = p_sites(i)%id
        exit
      end if
    end do
    if (site_index == -1) then
      print *, "Error: Position not found in lattice."
      print *, "Position:", rpos
      stop
    end if
    re = site_index
  end function get_site_index_from_rpos

  subroutine get_dpos_from_uc_index(ind,d,length,dpos)
    implicit none
    ! This function returns the integer position from the unit cell index
    ! d is the dimension of the position vector, length is the length of the lattice in each dimension
    ! dpos is the position vector to be returned, where
    ! each component tops [L1, L2, L3, ...](1-based) or [L1-1, L2-1, L3-1, ...](0-based)
    ! ind is the index of the unit cell in the lattice
    ! ind = 1 + (dpos(1)) + L1*((dpos(2)) + L2*((dpos(3)) + ...))) (0-based)
    ! ind = 1 + (dpos(1)-1) + L1*((dpos(2)-1) + L2*((dpos(3)-1) + ...))) (1-based)
    ! this subroutine is self-contained, it does not depend on the lattice type or outside variables
    
    integer, intent(in) :: d, ind, length(d)
    integer,intent(inout) :: dpos(d)
    integer :: i,ind_temp
    
    ! Initialize dpos to zero
    dpos = 0
    ind_temp = ind - 1 ! Convert to zero-based index
    do i = 1, d
      !dpos(i) = mod(ind_temp, length(i)) + 1 ! +1 to convert to 1-based index
      dpos(i) = mod(ind_temp, length(i)) ! No +1 here, as dpos is 0-based
      !ind_temp = (ind_temp - dpos(i)+1) / length(i) ! Update ind_temp for the next dimension(1-based)
      ind_temp = (ind_temp - dpos(i)) / length(i) ! Update ind_temp for the next dimension(0-based)
    end do
    
  end subroutine get_dpos_from_uc_index

  subroutine get_uc_index_from_dpos(dpos,ind)
    implicit none
    ! This function returns the unit cell index from the integer position
    ! dpos is the position vector to be input, where
    ! each component tops [L1-1, L2-1, L3-1, ...](1-based)
    ! ind is the index of the unit cell in the lattice to be returned
    ! ind = 1 + (dpos(1)) + L1*((dpos(2)) + L2*((dpos(3)) + ...))) (0-based)
    ! ind = 1 + (dpos(1)-1) + L1*((dpos(2)-1) + L2*((dpos(3)-1) + ...))) (1-based)
    ! this subroutine is self-contained, it does not depend on the lattice type or outside variables
    
    integer, intent(in) :: dpos(:)
    integer, intent(out) :: ind
    integer :: dpos_copy(size(dpos))
    integer :: i
    dpos_copy = dpos
    
    do i = 1,size(dpos)
      if(dpos_copy(i) < 0) dpos_copy(i) = dpos_copy(i) + lat%dlength(i) ! Handle negative positions
      dpos_copy(i) = mod(dpos_copy(i), lat%dlength(i)) ! Ensure dpos is within bounds (0-based)
    end do
    ind = 1 ! Start from 1 for 1-based index
    do i = size(dpos), 1, -1
      ind = ind + dpos_copy(i) * product(lat%dlength(1:i-1), mask=(i>1))
    end do
    
  end subroutine get_uc_index_from_dpos

  subroutine get_relative_index(ind, dpos1, dpos2)
    implicit none
    ! This function returns the relative index from dpos1 to dpos2
    ! considering periodic boundary conditions
    ! dpos1 and dpos2 are the position vectors to be input, where
    ! each component tops [L1-1, L2-1, L3-1, ...](0-based)
    ! ind is the relative index to be returned
    ! ind = 1 + (dpos_rel(1)) + L1*((dpos_rel(2)) + L2*((dpos_rel(3)) + ...))) (0-based)
    ! ind = 1 + (dpos_rel(1)-1) + L1*((dpos_rel(2)-1) + L2*((dpos_rel(3)-1) + ...))) (1-based)
    ! this subroutine is self-contained, it does not depend on the lattice type or outside variables

    integer, intent(in) :: dpos1(lat%dim), dpos2(lat%dim)
    integer, intent(out) :: ind
    integer :: i
    integer :: dpos_rel(size(dpos1))
    
      ! Calculate relative position considering periodic boundary conditions
      do i = 1, size(dpos1)
        dpos_rel(i) = mod(dpos2(i) - dpos1(i) + lat%dlength(i), lat%dlength(i))
      end do
    
    ! Now calculate the index from the relative position
    call get_uc_index_from_dpos(dpos_rel, ind)

  end subroutine get_relative_index
end module lattice
