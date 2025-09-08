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
  implicit none
  integer ::  n_obs,total_width,n_bins,n_meas_per_bin,field_width,n_cores,MPI_nblock,MPI_one_block
  integer :: lat_dim,i_block
  integer,allocatable :: lat_dlength(:)
  character(len = 20) :: output_format
  character(len=100) :: prefix1 ='data/data/', prefix2='data/out_files/', &
  & nml_file = 'reading_guide.nml',data_file = 'out_core'
  type :: observable
    character(len=40) :: name
    integer :: kind,width,off_lo,off_hi
    real(8), allocatable :: block_mean(:,:,:,:) ! (n_bins,mpi_one_block,obs%width,MPI_nblock)
    real(8), allocatable :: out_mean(:,:) ! (width,mpi_nblock)
    real(8), allocatable :: out_err(:,:)  ! (width,mpi_nblock)
  end type observable
  type(observable), allocatable :: obs_list(:)
  namelist /basics/ n_obs,total_width,n_bins,n_meas_per_bin,field_width,n_cores,MPI_nblock,MPI_one_block,output_format,lat_dim
  namelist /Lat/ lat_dlength

  call readin_guide_and_init()

  call readin_block_mean()

  call compute_out_mean_and_err()

  do i_block = 1,MPI_nblock
    call output_results(i_block)
  end do
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
    integer :: i, lun, ios,j,i_block,i_core,line,i_obs,myid
    character(len=100) :: ci_myid,temp_file,ci_bin
    write(ci_bin,'(I4)') n_bins

    do myid=0,n_cores-1
      i_block = myid / MPI_one_block + 1
      i_core = mod(myid,MPI_one_block) + 1
      write (ci_myid, '(i4)') myid
      temp_file = trim(prefix1)//trim(data_file)//trim(adjustl(ci_myid))//'.dat'
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
!------------------------------------------------------!
!--------------------Statistics------------------------!
!------------------------------------------------------!
  subroutine compute_out_mean_and_err()
    implicit none
    integer :: i_obs,i_core,i_block,i_bin
    real(8), allocatable :: temp_sum(:)

    do i_obs = 1,n_obs

      allocate(temp_sum(obs_list(i_obs)%width))
      ! compute the i_block mean: sum(data)/(n_bins*MPI_one_block)
      do i_block = 1,MPI_nblock
        temp_sum =0d0
        do i_core = 1 , MPI_one_block
          temp_sum = temp_sum + sum(obs_list(i_obs)%block_mean(:,i_core,:,i_block),1)
        end do
        obs_list(i_obs)%out_mean(:,i_block) = temp_sum/real(MPI_one_block*n_bins)
      end do
      ! compute the error bar: sqrt(sum((data-mean)^2)/(n_bins*MPI_one_block-1)/(n_bins*MPI_one_block))
      do i_block = 1,MPI_nblock
        temp_sum =0d0
        do i_core = 1 , MPI_one_block
          do i_bin = 1,n_bins
            temp_sum = temp_sum + &
            & sum((obs_list(i_obs)%block_mean(i_bin,i_core,:,i_block) - obs_list(i_obs)%out_mean(:,i_block))**2,1)
          end do
        end do
        obs_list(i_obs)%out_err(:,i_block) = sqrt(temp_sum/real(MPI_one_block*n_bins-1)/real(MPI_one_block*n_bins))
      end do
      deallocate(temp_sum)
    end do
  end subroutine

!------------------------------------------------------!
!-----------------------Output-------------------------!
!------------------------------------------------------!

  subroutine output_results(i_block)
    integer,intent(in) :: i_block
    integer :: i_obs,line,i
    character(len=100) ::ci1,filename
    namelist /file/ filename
    write(ci1,'(I4)') i_block
    ci1 = trim(prefix2)//'b'//trim(adjustl(ci1))//'.nml'
    open(unit=30,file=ci1,status='old',action='read')
    read(30,nml=file)
    close(30)
    open(unit=30,file=trim(filename),status='old',action='write',position='append')
    write(30,'(A15,A6,A18,A18)') 'Name           ','Index','Average','Error'
    write(30,'(A)') ' '
    do i_obs=1,n_obs
      if(obs_list(i_obs)%width>1) then
        write(30,'(A)') '#-----------------------------------------------------------'
      end if
      do line = 1,obs_list(i_obs)%width
        write(30,'(A15,1I6,1'//trim(output_format)//',1'//trim(output_format)//')') adjustl(obs_list(i_obs)%name),line,&
        & obs_list(i_obs)%out_mean(line,i_block),obs_list(i_obs)%out_err(line,i_block)
      end do
      write(30,*)' '
    end do
    close(30)
  end subroutine output_results

end program output_new