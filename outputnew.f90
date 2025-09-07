program output_new
  integer ::  n_obs,total_width,n_bins,n_meas_per_bin,field_width,n_cores,MPI_nblock,MPI_one_block
  character(len = 20) :: output_format
  character(len=100) :: prefix ='data/data/', nml_file = 'reading_guide.nml',data_file = 'out_core'
  type :: observable
    character(len=40) :: name
    integer :: kind,width,off_lo,off_hi
    real(8), allocatable :: block_mean(:,:,:,:) ! (n_bins,mpi_one_block,obs%width,MPI_nblock)
    real(8), allocatable :: out_mean(:,:) ! (width,mpi_nblock)
    real(8), allocatable :: out_err(:,:)  ! (width,mpi_nblock)
  end type observable
  type(observable), allocatable :: obs_list(:)
  namelist /basics/ n_obs,total_width,n_bins,n_meas_per_bin,field_width,n_cores,MPI_nblock,MPI_one_block,output_format

  call readin_guide_and_init()

  call readin_block_mean()

  call compute_out_mean_and_err()

contains
  subroutine readin_guide_and_init()
    integer :: i, lun, ios
    character(len = 40) :: name
    integer :: kind,width,offset_lo,offset_hi
    namelist /obs/ name,kind,width,offset_lo,offset_hi
    open(unit=10,file=trim(prefix)//trim(nml_file),status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'cannot open reading guide file'
    read(10,nml=basics,iostat=ios)
    if (ios /= 0) stop 'error reading reading guide file'
    close(10)
    allocate(obs_list(n_obs))
    open(unit=10,file=trim(prefix)//trim(nml_file),status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'cannot open obs_list file'
    rewind(10)
    do i=1,n_obs
      read(10,nml=obs,iostat=ios)
      if (ios /= 0) stop 'error reading obs_list file'
      obs_list(i)%name = name
      obs_list(i)%kind = kind
      obs_list(i)%width = width
      obs_list(i)%off_lo = offset_lo
      obs_list(i)%off_hi = offset_hi
      allocate(obs_list(i)%block_mean(n_bins,MPI_one_block,obs_list(i)%width,MPI_nblock))
      allocate(obs_list(i)%out_mean(obs_list(i)%width,MPI_nblock))
      allocate(obs_list(i)%out_err(obs_list(i)%width,MPI_nblock))
      obs_list(i)%block_mean = 0.0d0
      obs_list(i)%out_mean = 0.0d0
      obs_list(i)%out_err = 0.0d0
    end do
  end subroutine readin_guide_and_init

  subroutine readin_block_mean()
    integer :: i, lun, ios,j,i_block,i_core,line,i_obs
    character(len=100) :: ci_myid,temp_file,ci_bin
    write(ci_bin,'(I4)') n_bins

    do myid=0,n_cores-1
      i_block = myid / MPI_one_block + 1
      i_core = mod(myid,MPI_one_block) + 1
      write (ci_myid, '(i4)') myid
      temp_file = trim(prefix)//trim(data_file)//trim(adjustl(ci_myid))//'.dat'
      open(unit=20,file=temp_file,status='old',action='read')
      !if (ios /= 0) stop 'cannot open data file'
      do i_obs=1,n_obs
        do line = 1,obs_list(i_obs)%width
          read(unit=20,FMT = '('//trim(adjustl(ci_bin))//trim(output_format)//')',iostat=ios) &
          & obs_list(i_obs)%block_mean(:,i_core,line,i_block)
        end do
      end do
    end do
    close(20)
  end subroutine readin_block_mean

  subroutine compute_out_mean_and_err()
    integer :: i_obs,i_core,i_block
    real(8), allocatable :: temp_sum(:)
    implicit none
    do i_obs = 1,n_obs
      allocate(temp_sum(obs_list(i_obs)%width))
      do i_block = 1,MPI_nblock
        temp_sum =0d0
        do i_core = 1 , MPI_one_block
          temp_sum = temp_sum + sum(obs_list(i_obs)%block_mean(:,i_core,:,i_block),1)
        end do
        obs_list(i_obs)%out_mean(:,i_block) = temp_sum/real(MPI_one_block*n_bins)
      end do

      deallocate(temp_sum)
    end do
  end subroutine
end program output_new