program output_new
!------------------------------------------------------!
  !this program readin the i_block mean data from all cores
  !and compute the final mean and error bar
  !then output the results into files
  ! prefix1: the address of bin average data, prefix2: the address of target, output folder
  ! nml_file: the name of reading guide file, data_file: the prefix of bin average data file
  ! all the information of main QMC code is passed through the reading guide file
  ! to output, we need to know the name of the target file, which is in 'b*.nml' file
  ! the target file should be opened in old mode, and the data will be appended to the end of the file
  !------------------------------------------------------!
  use nrtype
  implicit none
  integer ::  n_obs,total_width,n_bins,n_meas_per_bin,field_width,n_cores,MPI_nblock,MPI_one_block
  integer :: lat_dim,i_block,size_of_dp
  integer,allocatable :: lat_dlength(:)
  character(len = 20) :: output_format
  character(len=100) :: prefix1 ='data/data/', prefix2='data/out_files/', &
  & nml_file = 'reading_guide.nml',data_file = 'out_core'
  real(dp),allocatable :: binary_data(:,:) ! (total_width,n_bins)

  type :: observable
    character(len=40) :: name
    integer :: kind,width,off_lo,off_hi
    real(dp), allocatable :: bin_data(:,:,:,:) !> (obs%width,n_bins,mpi_one_block=n_cores,MPI_nblock)
    real(dp), allocatable :: core_mean(:,:,:) !>(obs%width,mpi_one_block=n_cores,MPI_nblock)
    real(dp), allocatable :: out_mean(:,:) !> (width,mpi_nblock)
    real(dp), allocatable :: out_err(:,:)  !> (width,mpi_nblock),assuming all bins are independent
    real(dp), allocatable :: large_err(:,:) !> (width,mpi_nblock), asuming bins in one block are correlated
  end type observable
  type(observable), allocatable :: obs_list(:)
  namelist /basics/ n_obs,total_width,n_bins,n_meas_per_bin,& 
  & field_width,n_cores,MPI_nblock,MPI_one_block,output_format,lat_dim,size_of_dp
  namelist /Lat/ lat_dlength

  call readin_guide_and_init()

  call readin_bin_data()

  call compute_out_mean_and_err()
  
  do i_block = 1,MPI_nblock
    call output_results(i_block)
  end do
  if(.true.) then
    print*, obs_list(1)%bin_data(:,:,1,1)
    print*, obs_list(1)%bin_data(:,:,2,1)
    print*, obs_list(1)%core_mean(:,:,1)
    print*, obs_list(1)%out_mean(:,1)
    print*, obs_list(1)%out_err(:,1)
    print*, obs_list(1)%large_err(:,1)
  end if
contains

!------------------------------------------------------!
!------------------initialization----------------------!
!------------------------------------------------------!
  subroutine readin_guide_and_init()
    integer :: i, lun, ios
    character(len = 40) :: name
    integer :: kind,width,offset_lo,offset_hi
    namelist /obs/ name,kind,width,offset_lo,offset_hi
    open(unit=10,file=trim(prefix1)//trim(nml_file),status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'cannot open reading guide file'
    read(10,nml=basics,iostat=ios)
    if (ios /= 0) stop 'error reading reading guide file'
    close(10)
    allocate(obs_list(n_obs))
    open(unit=10,file=trim(prefix1)//trim(nml_file),status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'cannot open obs_list file'
    rewind(10)
    allocate(binary_data(total_width,n_bins))
    do i=1,n_obs
      read(10,nml=obs,iostat=ios)
      if (ios /= 0) stop 'error reading obs_list file'
      obs_list(i)%name = name
      obs_list(i)%kind = kind
      obs_list(i)%width = width
      obs_list(i)%off_lo = offset_lo
      obs_list(i)%off_hi = offset_hi
      allocate(obs_list(i)%bin_data(obs_list(i)%width,n_bins,MPI_one_block,MPI_nblock))
      allocate(obs_list(i)%core_mean(obs_list(i)%width,MPI_one_block,MPI_nblock))
      allocate(obs_list(i)%out_mean(obs_list(i)%width,MPI_nblock))
      allocate(obs_list(i)%out_err(obs_list(i)%width,MPI_nblock))
      allocate(obs_list(i)%large_err(obs_list(i)%width,MPI_nblock))
      obs_list(i)%bin_data = 0.0d0
      obs_list(i)%core_mean = 0.0d0
      obs_list(i)%out_mean = 0.0d0
      obs_list(i)%out_err = 0.0d0
      obs_list(i)%large_err = 0.0d0
    end do
  end subroutine readin_guide_and_init

  subroutine readin_bin_data()
    integer ::  ios,i_core,line,i_obs,myid,i_bin
    character(len=100) :: ci_myid,temp_file,ci_bin
    write(ci_bin,'(I4)') n_bins
    if(.not.allocated(binary_data)) then
      allocate(binary_data(total_width,n_bins))
    end if
    ! binary_data(total_width,n_bins) read in the out_core*.bin file from each core
    do myid=0,n_cores-1
      ! readin from binary file
      binary_data = 0.0d0
      i_block = myid / MPI_one_block + 1
      i_core = mod(myid,MPI_one_block) + 1
      write (ci_myid, '(i4)') myid
      temp_file = trim(prefix1)//trim(data_file)//trim(adjustl(ci_myid))//'.bin'
      
      open(unit=20,file=temp_file,status='old',access='stream', form='unformatted')
      if (ios /= 0) stop 'cannot open binary data file'
      read(20) binary_data
      close(20)
      ! distribute data to obs_list
      !do i_obs=1,n_obs
      !  do line = 1,obs_list(i_obs)%width
      !    read(unit=20,FMT = '('//trim(adjustl(ci_bin))//trim(output_format)//')',iostat=ios) &
      !    & obs_list(i_obs)%bin_data(:,i_core,line,i_block)
      !  end do 
      !end do
      do i_obs=1,n_obs
        do i_bin = 1,n_bins
          obs_list(i_obs)%bin_data(:,i_bin,i_core,i_block) = &
          & binary_data(obs_list(i_obs)%off_lo:obs_list(i_obs)%off_hi,i_bin)
        end do
      end do
    end do
    
  end subroutine readin_bin_data
!------------------------------------------------------!
!--------------------Statistics------------------------!
!------------------------------------------------------!
  subroutine compute_out_mean_and_err()
    implicit none
    integer :: i_obs,i_core,i_block,i_bin
    real(dp), allocatable :: temp_sum(:,:) ! > (obs%width, mpi_one_block=n_cores)

    do i_obs = 1,n_obs
      ! compute the core mean: sum(data)/n_bins
      obs_list(i_obs)%core_mean = sum(obs_list(i_obs)%bin_data,dim=2)/real(n_bins)
      allocate(temp_sum(obs_list(i_obs)%width,MPI_one_block))
      ! compute the i_block mean: sum(data)/(n_bins*MPI_one_block) and core mean: sum(data)/n_bins
      do i_block = 1,MPI_nblock
        temp_sum =0d0
        do i_core = 1 , MPI_one_block
          do i_bin = 1,n_bins
            temp_sum(:,i_core) = temp_sum(:,i_core) + obs_list(i_obs)%bin_data(:,i_bin,i_core,i_block)
          end do
        end do
        obs_list(i_obs)%out_mean(:,i_block) = sum(temp_sum,dim=2)/real(MPI_one_block*n_bins)
      end do

      ! compute the error bar: sqrt(sum((data-mean)^2)/(n_bins*MPI_one_block-1)/(n_bins*MPI_one_block))
      do i_block = 1,MPI_nblock
        temp_sum =0d0
        do i_core = 1 , MPI_one_block
          do i_bin = 1,n_bins
            temp_sum(:,i_core) = temp_sum(:,i_core) + &
            & (obs_list(i_obs)%bin_data(:,i_bin,i_core,i_block) - obs_list(i_obs)%out_mean(:,i_block))**2

          end do
        end do
        obs_list(i_obs)%out_err(:,i_block) = sqrt(sum(temp_sum,dim=2)/real(MPI_one_block*n_bins-1)/real(MPI_one_block*n_bins))
      end do

      ! compute the large error bar: sqrt(sum((core_mean-mean)^2)/(MPI_one_block-1)/(MPI_one_block))
      do i_block = 1,MPI_nblock
        temp_sum =0d0
        do i_core = 1 , MPI_one_block
          temp_sum(:,i_core) = temp_sum(:,i_core) + & 
          & (obs_list(i_obs)%core_mean(:,i_core,i_block) - obs_list(i_obs)%out_mean(:,i_block))**2
        end do
        obs_list(i_obs)%large_err(:,i_block) = sqrt(sum(temp_sum,dim=2)/real(MPI_one_block-1)/real(MPI_one_block))
      end do
      deallocate(temp_sum)
    end do
  end subroutine

!------------------------------------------------------!
!-----------------------Output-------------------------!
!------------------------------------------------------!

  subroutine output_results(i_block)
    integer,intent(in) :: i_block
    integer :: i_obs,line,i,count
    character(len=100) ::ci1,filename
    character(len = 100):: readin_line
    logical :: found = .false.
    namelist /file/ filename
    write(ci1,'(I4)') i_block
    ci1 = trim(prefix2)//'b'//trim(adjustl(ci1))//'.nml'
    open(unit=30,file=ci1,status='old',action='read')
    read(30,nml=file)
    close(30)
    !--------------open the file and truncate it------------------!
    open(unit=30,file=trim(filename),status='old',action='readwrite', &
       form='formatted', access='sequential')
    rewind(30)
    count = 0
    do while (.not.found)
      count = count + 1
      read(30,'(A)') readin_line
      print*,trim(adjustl(readin_line))
      if(trim(adjustl(readin_line)) == '#--observables--#') then
        found = .true.
        endfile(30)    
        exit
      end if
    end do
    close(30)
    !--------------append the data to the end of the file---------!

    open(unit=30,file=trim(filename),status='old',action='readwrite',position='append')
    write(30,*)' '
    write(30,'(A15,A6,A18,A18,A18)') 'Name           ','Index','Average','Error', 'LError'
    write(30,'(A)') ' '
    do i_obs=1,n_obs
      if(obs_list(i_obs)%width>1) then
        write(30,'(A)') '#-----------------------------------------------------------'
      end if
      do line = 1,obs_list(i_obs)%width
        write(30,'(A15,1I6,1'//trim(output_format)//',1'//trim(output_format)//',1'//trim(output_format)//')') & 
        & adjustl(obs_list(i_obs)%name),line,&
        & obs_list(i_obs)%out_mean(line,i_block),obs_list(i_obs)%out_err(line,i_block), obs_list(i_obs)%large_err(line,i_block)
        if(obs_list(i_obs)%width > 1) then
          !write(*,'(1'//trim(output_format)//')') obs_list(i_obs)%out_err(line,i_block)
        end if
      end do
      write(30,*)' '
    end do
    close(30)
  end subroutine output_results

end program output_new