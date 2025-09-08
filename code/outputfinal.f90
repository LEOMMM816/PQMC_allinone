program outputfinal
   use input
   implicit none
   logical :: not_correlated = .true.
   character*100 :: ci, ci1, ci2, ci3 ! name of file
   character*100 :: str , jobname
   character(len=40) :: path = 'data/data/'
   integer :: i_block,obs
   REAL(8), PARAMETER :: PI=3.141592653589793238462643383279502884197d0
   real(8),allocatable :: obs_block_aver(:, :), obs_block_err(:,:), corf_block_aver(:,:,:), corf_block_err(:,:,:)
   real(8), allocatable :: corfk_block_aver(:,:,:), corfk_block_err(:,:,:)
   real(8),allocatable :: ctau_block_aver(:,:,:), ctau_block_err(:,:,:),cw_block_aver(:,:,:), cw_block_err(:,:,:)
   !character*30, allocatable :: obs_name(:)
   call get_command_argument(1, jobname)
   call init(0)
   if(.not.allocated(obs_block_aver)) allocate (obs_block_aver(MPI_nblock, n_obs), obs_block_err(MPI_nblock, n_obs))
   if(.not.allocated(corf_block_aver)) allocate (corf_block_aver(corf_length,MPI_nblock, n_corf), & 
   & corf_block_err(corf_length,MPI_nblock, n_corf))
   if(.not.allocated(corfk_block_aver)) allocate (corfk_block_aver(corf_length,MPI_nblock, n_corf), &
   & corfk_block_err(corf_length,MPI_nblock, n_corf))
   if(.not.allocated(ctau_block_aver)) allocate (ctau_block_aver(ntau+1,MPI_nblock, n_ctau), &
   & ctau_block_err(ntau+1,MPI_nblock, n_ctau))
   if(.not.allocated(cw_block_aver)) allocate (cw_block_aver(ntau+1,MPI_nblock, n_ctau), &
   & cw_block_err(ntau+1,MPI_nblock, n_ctau))
   call MPI_output_final()
   call MPI_output_final_corf()
   do i_block = 1, MPI_nblock
      call init_MPI(i_block-1)
      call mpi_output_inforsheet(i_block)
   end do
CONTAINS
   subroutine MPI_output_final()

      real(8) :: obs_block(nbin_per_core*MPI_one_block, MPI_nblock, n_obs)
      real(8), allocatable :: obs_core_aver(:) ! (n_core)
      integer :: i_core
      integer :: total_line
      do obs = 1, n_obs
         do i_block = 1, MPI_nblock

            ci = obs_name(obs)
#ifdef MPI
            write (ci2, '(1i4)') i_block - 1
            str = 'o_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'.dat'
#else
            str = 'o_'//trim(adjustl(ci))//'.dat'
#endif
            str = trim(adjustl(path))//trim(adjustl(str))

            open (80 + obs, file=str, status='old', position='REWIND')
            total_line = GetFileN(80 + obs)
            close (80 + obs)
            if (.not. allocated(obs_core_aver)) allocate (obs_core_aver(total_line))

            open (80 + obs, file=str, status='old', position='REWIND')
            WRITE (ci, '(1i4)') nbin_per_core

            do i_core = 0, total_line - 1
      read (80 + obs, '('//trim(adjustl(ci))//'f18.12)') & 
      & obs_block(nbin_per_core*i_core + 1:nbin_per_core*(i_core + 1), i_block, obs)
         obs_core_aver(i_core + 1) = sum(obs_block(nbin_per_core*i_core + 1:nbin_per_core*(i_core + 1), i_block, obs))/nbin_per_core
            end do
            close (80 + obs)
            !obs_block(:,i_block,8) = abs(obs_block(:,i_block,8) - obs_block(:,i_block,4))
            obs_block_aver(i_block, obs) = sum(obs_block(1:nbin_per_core * total_line, i_block, obs))/(nbin_per_core*total_line)

            if(not_correlated) then
            obs_block_err(i_block,obs) = sum((obs_block(1:nbin_per_core * total_line,i_block,obs)-obs_block_aver(i_block,obs))**2)
            obs_block_err(i_block,obs) = sqrt(obs_block_err(i_block, obs))/sqrt(real(nbin_per_core*total_line)&
            & * (nbin_per_core*total_line-1))
            else
            ! if bins in one core are correlated
            obs_block_err(i_block, obs) = sum((obs_core_aver(:) - obs_block_aver(i_block, obs))**2) 
            obs_block_err(i_block, obs) = sqrt(obs_block_err(i_block, obs))/sqrt(real(total_line)*real(total_line-1))
            end if
            deallocate (obs_core_aver)
         end do
      end do
      if(.false.) then
      open (100, file=trim(adjustl(jobname))//'aver.dat', status='replace')
      write (ci, '(1i4)') n_obs
      write (100, '('//trim(adjustl(ci))//'A18)') obs_name
      do i_block = 1, MPI_nblock
         write (100, '('//trim(adjustl(ci))//'f18.12)') obs_block_aver(i_block, :)
      end do
      close (100)

      open (100, file=trim(adjustl(jobname))//'err.dat', status='replace')
      write (100, '('//trim(adjustl(ci))//'A18)') obs_name
      write (ci, '(1i4)') n_obs
      do i_block = 1, MPI_nblock
         write (100, '('//trim(adjustl(ci))//'f18.12)') obs_block_err(i_block, :)
      end do
      close (100)
      end if
   end subroutine


   subroutine MPI_output_final_corf()

      real(8) :: corf_block(corf_length,nbin_per_core*MPI_one_block, MPI_nblock, n_corf)
      real(8) :: corfk_block(corf_length,nbin_per_core*MPI_one_block, MPI_nblock, n_corf)
      real(8) :: ctau_block(ntau+1,nbin_per_core*MPI_one_block, MPI_nblock, n_ctau)
      real(8) :: cw_block(ntau+1,nbin_per_core*MPI_one_block, MPI_nblock, n_ctau)
      integer :: fucked_up_list(1) = (/0/)
      integer :: i_core, i_block,corf,i_site,s1,s2,bin,dbin
      integer :: total_line
      real(8) :: temp
      complex(8), allocatable :: corfk_phase_inv(:,:)
         ! read corfk_data from file
      do corf = 1, n_corf
         ci = corf_name(corf)
      do i_block = 1, MPI_nblock
         do i_core = 1, MPI_one_block
#ifdef MPI
            write (ci1, '(1i4)') i_core - 1
            write (ci2, '(1i4)') i_block - 1
            str = 'ck_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'c'//trim(adjustl(ci1))//'.dat'
#else
            str = 'ck_'//trim(adjustl(ci))//'.dat'
#endif
            str = trim(adjustl(path))//trim(adjustl(str))
            open (80 + corf, file=str, status='old', position='REWIND')
            do i_site = 1, corf_length
               WRITE (ci3, '(1i4)') nbin_per_core
               read (80 + corf, '('//trim(adjustl(ci3))//'f18.12)') & 
               & corfk_block(i_site,(i_core-1) * nbin_per_core + 1 : (i_core) * nbin_per_core, i_block,corf)
            end do
            close (80 + corf)
                  
         end do

         do i_site = 1, corf_length
         corfk_block_aver(i_site,i_block,corf) = sum(corfk_block(i_site,:, i_block,corf))/(nbin_per_core*MPI_one_block)
         corfk_block_err(i_site,i_block,corf) = sum((corfk_block(i_site,:,i_block,corf)-corfk_block_aver(i_site,i_block,corf))**2)
         corfk_block_err(i_site,i_block,corf) = sqrt(corfk_block_err(i_site,i_block, corf))/sqrt(real(nbin_per_core*MPI_one_block)&
         & * (nbin_per_core*MPI_one_block-1))
         
         end do

      end do
   end do
   ! read corf_data from file or Fourier transform from corfk_data
   if(.not. allocated(corfk_phase_inv)) allocate (corfk_phase_inv(corf_length, corf_length))
   do s1 = 1, corf_length
      ! s1  is the index of k , - pi to (corf_length - 1)/La * pi
         do s2 = 1, corf_length
            ! s2 is the index of r, 0 to (corf_length - 1)
            corfk_phase_inv(s2, s1) = exp(complex(0d0, -2d0 * pi * (s1 - 1-La/2) * (s2 - 1) / La))/corf_length
         end do
   end do

      do corf = 1, n_corf
         ci = corf_name(corf)
         
         do i_block = 1, MPI_nblock
            if(findloc(fucked_up_list,corf,dim = 1) == 0) then
            do i_core = 1, MPI_one_block
#ifdef MPI
               write (ci1, '(1i4)') i_core - 1
               write (ci2, '(1i4)') i_block - 1
               str = 'c_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'c'//trim(adjustl(ci1))//'.dat'
#else
               str = 'c_'//trim(adjustl(ci))//'.dat'
#endif
               str = trim(adjustl(path))//trim(adjustl(str))
               open (80 + corf, file=str, status='old', position='REWIND')
               do i_site = 1, corf_length
                  WRITE (ci3, '(1i4)') nbin_per_core
                  read (80 + corf, '('//trim(adjustl(ci3))//'f18.12)') & 
                  & corf_block(i_site,(i_core-1) * nbin_per_core + 1 : (i_core) * nbin_per_core, i_block,corf)
               end do
               close (80 + corf)
                     
            end do
         else 
               corf_block(:,:,i_block, corf) = real(matmul(corfk_phase_inv,corfk_block(:,:, i_block,corf)))
         end if
         do i_site = 1, corf_length
            corf_block_aver(i_site,i_block,corf) = sum(corf_block(i_site,:, i_block,corf))/(nbin_per_core*MPI_one_block)
            corf_block_err(i_site,i_block,corf) = sum((corf_block(i_site,:,i_block,corf)-corf_block_aver(i_site,i_block,corf))**2)
            corf_block_err(i_site,i_block,corf) = sqrt(corf_block_err(i_site,i_block, corf))/sqrt(real(nbin_per_core*MPI_one_block)&
            & * (nbin_per_core*MPI_one_block-1))
            
            end do
      end do ! end for i_block
      end do ! end for corf

      ! read ctau_data from file
      do corf = 1, n_ctau
         ci = ctau_name(corf)
      do i_block = 1, MPI_nblock
         do i_core = 1, MPI_one_block
#ifdef MPI
            write (ci1, '(1i4)') i_core - 1
            write (ci2, '(1i4)') i_block - 1
            str = 'ct_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'c'//trim(adjustl(ci1))//'.dat'
#else
            str = 'ct_'//trim(adjustl(ci))//'.dat'
#endif
            str = trim(adjustl(path))//trim(adjustl(str))
            open (80 + corf, file=str, status='old', position='REWIND')
            do i_site = 1, ntau+1
               WRITE (ci3, '(1i4)') nbin_per_core
               read (80 + corf, '('//trim(adjustl(ci3))//'f18.12)') & 
               & ctau_block(i_site,(i_core-1) * nbin_per_core + 1 : (i_core) * nbin_per_core, i_block,corf)
            end do
            close (80 + corf)
         
         end do
      end do
      end do

      ! read cw_data from file
      do corf = 1, n_ctau
         ci = ctau_name(corf) 
      do i_block = 1, MPI_nblock
         do i_core = 1, MPI_one_block
#ifdef MPI
            write (ci1, '(1i4)') i_core - 1
            write (ci2, '(1i4)') i_block - 1
            str = 'cw_'//trim(adjustl(ci))//'_b'//trim(adjustl(ci2))//'c'//trim(adjustl(ci1))//'.dat'
#else
            str = 'cw_'//trim(adjustl(ci))//'.dat'
#endif
            str = trim(adjustl(path))//trim(adjustl(str))
            open (80 + corf, file=str, status='old', position='REWIND')
            do i_site = 1, ntau+1
               WRITE (ci3, '(1i4)') nbin_per_core
               read (80 + corf, '('//trim(adjustl(ci3))//'f18.12)') & 
               & cw_block(i_site,(i_core-1) * nbin_per_core + 1 : (i_core) * nbin_per_core, i_block,corf)
            end do
            close (80 + corf)
         end do
         end do
      end do
      
      ! pick out ctau value > 1
      do corf = 1, n_ctau
         do i_block = 1, MPI_nblock
            do i_core = 1, MPI_one_block
               do bin = ((i_core-1) * nbin_per_core + 1) , (i_core) * nbin_per_core

                  temp = maxval((ctau_block(:,bin,i_block,corf)))
            if(temp > 1d0) then
               dbin = 1
               do while (temp > 1d0)
                  print*,'corf',corf,'block',i_block,'i_core',i_core,'temp',temp
                  ctau_block(:,bin,i_block,corf) = ctau_block(:,modi(bin-dbin,(i_core) * nbin_per_core),i_block,corf)
                  cw_block(:,bin,i_block,corf) = cw_block(:,modi(bin-dbin,(i_core) * nbin_per_core),i_block,corf)
                  temp = maxval((ctau_block(:,bin,i_block,corf)))
                  dbin = dbin + 1
               end do
            end if
         end do;end do;end do;end do
      
      ! calculate aver and err for ctau and cw
      do corf = 1 , n_ctau
      do i_block = 1, MPI_nblock
         do i_site = 1, ntau+1
         ctau_block_aver(i_site,i_block, corf) = sum(ctau_block(i_site,:, i_block,corf))/(nbin_per_core*MPI_one_block)
         ctau_block_err(i_site,i_block, corf) = sum((ctau_block(i_site,:,i_block,corf)-ctau_block_aver(i_site,i_block,corf))**2)
         ctau_block_err(i_site,i_block, corf) = sqrt(ctau_block_err(i_site,i_block, corf))/sqrt(real(nbin_per_core*MPI_one_block)&
            & * (nbin_per_core*MPI_one_block-1))

         cw_block_aver(i_site,i_block, corf) = sum(cw_block(i_site,:, i_block,corf))/(nbin_per_core*MPI_one_block)
         cw_block_err(i_site,i_block, corf) = sum((cw_block(i_site,:,i_block,corf)-cw_block_aver(i_site,i_block,corf))**2)
         cw_block_err(i_site,i_block, corf) = sqrt(cw_block_err(i_site,i_block, corf))/sqrt(real(nbin_per_core*MPI_one_block)&
            & * (nbin_per_core*MPI_one_block-1))

         end do
      end do
      end do
      ! write to file corfaver and corferr
      path = 'data/'
      ! begin write
      do corf = 1, n_corf
         ci = corf_name(corf)
         write(ci1, '(1i4)') MPI_nblock
         open (100, file=trim(adjustl(path))//'corfaver_'//trim(adjustl(ci))//'.dat', status='replace')
         do i_site = 1, corf_length
            write (100, '('//trim(adjustl(ci1))//'f18.6)') corf_block_aver(i_site, :, corf)
         end do
         close (100)
         open (100, file=trim(adjustl(path))//'corferr_'//trim(adjustl(ci))//'.dat', status='replace')
         do i_site = 1, corf_length
            write (100, '('//trim(adjustl(ci1))//'f18.6)') corf_block_err(i_site, :, corf)
         end do
         close (100)
      end do



   ! write corfk to file corfkaver and corfkerr
   do corf = 1, n_corf
      ci = corf_name(corf)

      write(ci1, '(1i4)') MPI_nblock
      open (100, file=trim(adjustl(path))//'corfkaver_'//trim(adjustl(ci))//'.dat', status='replace')
      do i_site = 1, corf_length
         write (100, '('//trim(adjustl(ci1))//'f18.6)') corfk_block_aver(i_site, :, corf)
      end do
      close (100)
      open (100, file=trim(adjustl(path))//'corfkerr_'//trim(adjustl(ci))//'.dat', status='replace')
      do i_site = 1, corf_length
         write (100, '('//trim(adjustl(ci1))//'f18.6)') corfk_block_err(i_site, :, corf)
      end do
      close (100)
   end do
        ! write ctau to file ctauaver and ctaueerr
   do corf = 1, n_ctau
      ci = ctau_name(corf)
      write(ci1, '(1i4)') MPI_nblock
      open (100, file=trim(adjustl(path))//'ctaver_'//trim(adjustl(ci))//'.dat', status='replace')
      do i_site = 1, ntau+1
         write (100, '('//trim(adjustl(ci1))//'f18.6)') ctau_block_aver(i_site, :, corf)
         
      end do
      close (100)
      open (100, file=trim(adjustl(path))//'cterr_'//trim(adjustl(ci))//'.dat', status='replace')
      do i_site = 1, ntau+1
         write (100, '('//trim(adjustl(ci1))//'f18.6)') ctau_block_err(i_site, :, corf)
      end do
      close (100)
   end do

   end subroutine
   Integer Function GetFileN(iFileUnit)
      Implicit None
      Integer, Intent(IN)::iFileUnit
      Integer::ios
      Character(Len=1)::cDummy
      GetFileN = 0
      Rewind (iFileUnit)
      Do
         Read (iFileUnit, *, ioStat=ioS) cDummy
         if (ioS /= 0) Exit
         GetFileN = GetFileN + 1
      End Do
   end function

   subroutine mpi_output_inforsheet(block)
      implicit none
      !! output all information of this simulation
      !! take the format of the La12Lb12Lt96dau... .out file in this folder
      integer :: block
      character(len=40) :: filename
      integer :: unit, i,i_corf,iterations
      write(ci, '(1i4)') block
      if(M > 0.0001d0) then
         ! phonon model
      write(filename, '(a,1i2,a,1i2,a,1i4,a,1f4.2,a,1f4.1,a,1i3,a)') 'Lx', La, 'Ly', Lb, 'Lt', ntime, 'dtau', delt, & 
      & 'lam', real(ep_parameter), 'Nup',nelec(1),'.out'
      else
         ! hubbard model
      write(filename, '(a,1i2,a,1i2,a,1i4,a,1f4.2,a,1f4.1,a,1i3,a)') 'Lx', La, 'Ly', Lb, 'Lt', ntime, 'dtau', delt, & 
      & 'U', U, 'Nup',nelec(1),'.out'
      end if
      call remove_spaces(filename,len(filename))
      iterations = (nbin_per_core*one_bin)*meas_interval + warmup
      open(newunit=unit, file=filename, status='replace', action='write')

  ! 写入文件头
      write(unit, '(a)') repeat('*', 48)
      write(unit, '(a)') repeat(' ', 1) // repeat('*', 47)
      write(unit, '(a)') repeat(' ', 11) // 'Holstein Creutz model' // repeat(' ', 11) // '*'
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
      write(unit, '(a, i15)') ' measurements in one_bin', one_bin
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
      write(unit, '(a, i15)') ' Lx    = ', La
      write(unit, '(a, i15)') ' Ly    = ', Lb
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
  do i = 1, n_obs
    write(unit, '(a15, f20.5, a, es20.6)') ' '//trim(obs_name(i)), obs_block_aver(block,i), &
                                         '      +/- ', obs_block_err(block,i)
  end do
  write(unit, '(a)') ''

  write(unit, '(a)') 'List of coorelation functions and errors'
      write(unit, '(a)') repeat(' ', 1) // 'Description: '
      write(unit, '(a)') repeat(' ', 3) // 'Most of the correlation functions are calcualted only on A sublattice'
      write(unit, '(a)') repeat(' ', 3) // 'The correlation functions with AB are calculated on across A and B sublattice'
      write(unit, '(a)') repeat(' ', 3) // 'given by <O_A(0)O_B(r)>'
  ! 写入实关联函数数据
  do i_corf = 1 , n_corf
      write(unit, '(a)') ' Average '//trim(corf_name(i_corf))//' correlation function'

      write(unit, '(a)') ''
      do i = 1, corf_length
            write(unit, '(a, 1i12, f20.8, es20.6)') ' '//trim(corf_name(i_corf))//'(r)', &
            i-1, corf_block_aver(i,block,i_corf), corf_block_err(i,block,i_corf)
      end do
      write(unit, '(a)') ''
   end do
   ! 写入k空间关联函数数据
   write(unit, '(a)') 'The numbering of k-space is from -pi to pi-\delta'
   write(unit, '(a)') ''
   do i_corf = 1 , n_corf
      if(corf_name(i_corf) == 'den-den') then
      write(unit, '(a)') ' Average k-space '//trim(corf_name(i_corf))//' correlation function = 1 - n(k)'
      else
         write(unit, '(a)') ' Average k-space '//trim(corf_name(i_corf))//' correlation function'
      end if
      write(unit, '(a)') ''
      do i = 1, corf_length
            write(unit, '(a, 1i12, f20.8, es20.6)') ' '//trim(corf_name(i_corf))//'(p)', &
            i-1-corf_length/2, corfk_block_aver(i,block,i_corf), corfk_block_err(i,block,i_corf)
      end do
      write(unit, '(a)') ''
   end do
   ! 写入时间关联函数数据
   write(unit, '(a)') 'The numbering of tau is from 0 to ntau'
   write(unit, '(a)') ''
   do i_corf = 1 , n_ctau
      write(unit, '(a)') ' Average time '//trim(ctau_name(i_corf))//' correlation function'
      write(unit, '(a)') ''
      do i = 1, ntau+1
            write(unit, '(a, 1i12, f20.8, es20.6)') ' '//trim(ctau_name(i_corf))//'(t)', &
            i-1, ctau_block_aver(i,block,i_corf), ctau_block_err(i,block,i_corf)
      end do
      write(unit, '(a)') ''
   end do
   ! 写入频率关联函数数据
   write(unit, '(a)') 'The numbering of w is from (-pi)+ to pi-'
   write(unit, '(a)') ''
   do i_corf = 1 , n_ctau
      write(unit, '(a)') ' Average frequency '//trim(ctau_name(i_corf))//' correlation function'
      write(unit, '(a)') ''
      do i = 1, ntau+1
            write(unit, '(a, 1i12, f20.8, es20.6)') ' '//trim(ctau_name(i_corf))//'(w)', &
            i-1-ntau/2, cw_block_aver(i,block,i_corf), cw_block_err(i,block,i_corf)
      end do
      write(unit, '(a)') ''
   end do
  ! 写入文件尾
  write(unit, '(a)') ''
  write(unit, '(a)') repeat('*', 48)
  write(unit, '(a)') repeat('*', 48)

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


integer function modi(n, p)
      integer, intent(in)::n, p

      modi = mod(n, p)
      if (n < 0) then
         modi = p + modi
         return

      else if (modi == 0) then
         modi = p
      end if
   end function
end program