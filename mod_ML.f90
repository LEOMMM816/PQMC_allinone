MODULE zxy_ML
   use ZXYP_SSH_server

   implicit none
! measurement parameter
   character*30 :: ci, ci2 ! name of file
   character*50 :: str
   integer :: ML_obs = 16, ML_sweep = 1, n_variable = 10, ml_group = 100
   integer :: ML_head = 1, debug_accept
   real(8) :: cumu_en
   real(8) :: J_trans, J_long, J_anti , J_ferro
! matrix declare
   real(8), allocatable :: ML_data(:, :) ! (n_sample,ML_obs)
   real(8), allocatable :: ML_Parameters(:, :), parameters(:)
! special type
   type :: bond_en
      integer :: bond
      real(8) :: en
   end type bond_en
! stack
   integer :: current_cluster, current_wl
   type(bond_en), allocatable :: wait_list(:)
   integer, allocatable :: cluster(:)
CONTAINS

   subroutine init_ml_parameter()
      implicit none
      integer :: i
      if (.not. allocated(parameters)) allocate (ML_parameters(n_variable, 16))
      if (.not. allocated(parameters)) allocate (parameters(n_variable))
      open (20, file='paraO05L8.dat', status='old', action='read')
      do i = 1, n_variable
         read (20, '(16E25.16)') ML_parameters(i, :)
         
      end do

      close (20)
      
      parameters = ML_parameters(:, MPI_block+1)
      print*,parameters
      J_anti = parameters(5)
      J_ferro = parameters(6)
      J_long = parameters(7)
      J_trans = parameters(8)

      !parameters = 0
   end subroutine
   subroutine ML_record()
      implicit none
      integer :: flag, i
      flag = 0
      if (.not. allocated(ML_data)) allocate (ML_data(ML_obs, ML_group))

      if (ml_head == 1) ML_data = 0d0

      call ML_cal_data(ML_data(:, ml_head))
      ml_head = ml_head + 1
      if (ml_head == ml_group + 1) then

         ml_head = 1
#ifdef MPI
         write (ci2, '(1i4)') MPI_block
         str = 'ML_b'//trim(adjustl(ci2))//'.dat'
#else
         str = 'ML.dat'
#endif
         str = 'ML_data/'//trim(adjustl(str))
         str = trim(adjustl(file_name))//trim(adjustl(str))
100      flag = flag + 1
         if (flag > 10) return
         inquire (file=str, exist=file_exist)
         inquire (file=str, exist=file_exist)
         if (.not. file_exist) then
            open (99, file=str, status='replace', ERR=100)
            close (99)
         end if
         open (99, file=str, status='old', position='append', ERR=100)
         !write(80+obs,'(f18.6)') beta(tem)
         WRITE (ci, '(1i4)') ML_obs
         do i = 1, ml_group
            write (99, '('//trim(adjustl(ci))//'f18.9)', ERR=101) ML_data(:, i)
101         cycle
         end do
         close (99)
      end if
   end subroutine

   subroutine ML_cal_data(ml_data_sub)
      implicit none
      real(8) :: ml_data_sub(ml_obs)
      integer :: site, site_1, site_2, bond_1, bond_2
      integer :: tau , i
      integer :: bond_near_1(8) , bond_near_2(8)
      ml_data_sub = 0d0
      !!data_1 : kenetic energy of ph
      do tau = 1, ntime
         ml_data_sub(1) = ml_data_sub(1) + sum((phonon_field(:, tau) - phonon_field(:, modi(tau + 1, ntime)))**2) &
         & /omega**2/delt**2/2d0
      END do
      ml_data_sub(1) = ml_data_sub(1)/Ns/ntime
      !!data_2 : potential energy of ph
      ml_data_sub(2) = sum(phonon_field**2)/2*D/Ns/ntime
      !!data_3: electron energy
      ml_data_sub(3) = ln_cw
      !!data_4: X_{i,l}^{2}
      ml_data_sub(4) = sum(phonon_field**2)/Ns/ntime
      !!data_5 : X_{i,l} * X_{i,l+1}

      do tau = 1, ntime
         ml_data_sub(5) = ml_data_sub(5) + sum(phonon_field(:, tau)*phonon_field(:, modi(tau + 1, ntime)))
      END do
      ml_data_sub(5) = ml_data_sub(5)/Ns/ntime

      !!data_6 : X_{i,l}^2 * X_{i,l+1}^2
      do tau = 1, ntime
         ml_data_sub(6) = ml_data_sub(6) + sum(phonon_field(:, tau)**2*phonon_field(:, modi(tau + 1, ntime))**2)
      END do
      ml_data_sub(6) = ml_data_sub(6)/Ns/ntime

      !!data_7 : X_{i,l}^{4}
      ml_data_sub(7) = sum(phonon_field**4)/Ns/ntime

      !!data_8 : X_{i,l}^{6}
      ml_data_sub(8) = sum(phonon_field**6)/Ns/ntime
      do site = 1, Ns
         call get_near_bond(site, bond_near_1)
         call get_near_bond(site + Ns, bond_near_2)
         !print*, 'site:', site
         !print*, 'bond_near_1:', bond_near_1
         !print*, 'bond_near_2:', bond_near_2
         !!data_9 : nearest anti interaction

         ml_data_sub(9) = ml_data_sub(9) + sum(phonon_field(site, :)*phonon_field(bond_near_1(1), :))
         ml_data_sub(9) = ml_data_sub(9) + sum(phonon_field(site, :)*phonon_field(bond_near_1(2), :))
         ml_data_sub(9) = ml_data_sub(9) + sum(phonon_field(site+Ns, :)*phonon_field(bond_near_2(1), :))
         ml_data_sub(9) = ml_data_sub(9) + sum(phonon_field(site+Ns, :)*phonon_field(bond_near_2(2), :))
         !!data_10 : nearest ferro interaction

         ml_data_sub(10) = ml_data_sub(10) + sum(phonon_field(site, :)*phonon_field(bond_near_1(3), :))
         ml_data_sub(10) = ml_data_sub(10) + sum(phonon_field(site, :)*phonon_field(bond_near_1(4), :))
         ml_data_sub(10) = ml_data_sub(10) + sum(phonon_field(site+Ns, :)*phonon_field(bond_near_2(3), :))
         ml_data_sub(10) = ml_data_sub(10) + sum(phonon_field(site+Ns, :)*phonon_field(bond_near_2(4), :))

         !!data_11 : X_{i,l} * X_{i+1,l} next nearest long

         !La direction
         ml_data_sub(11) = ml_data_sub(11) + sum((phonon_field(site, :)*phonon_field(bond_near_1(5), :)))
         ml_data_sub(11) = ml_data_sub(11) + sum((phonon_field(site, :)*phonon_field(bond_near_1(6), :)))
         !Lb direction
         ml_data_sub(11) = ml_data_sub(11) + sum((phonon_field(site+Ns, :)*phonon_field(bond_near_2(5), :)))
         ml_data_sub(11) = ml_data_sub(11) + sum((phonon_field(site+Ns, :)*phonon_field(bond_near_2(6), :)))
         !!data_12: X_{i,l} * X_{i+1,l} next nearest trans interaction

         !La direction
         ml_data_sub(12) = ml_data_sub(12) + sum((phonon_field(site, :)*phonon_field(bond_near_1(7), :)))
         ml_data_sub(12) = ml_data_sub(12) + sum((phonon_field(site, :)*phonon_field(bond_near_1(8), :)))

         !Lb direction
         ml_data_sub(12) = ml_data_sub(12) + sum((phonon_field(site+Ns, :)*phonon_field(bond_near_2(7), :)))
         ml_data_sub(12) = ml_data_sub(12) + sum((phonon_field(site+Ns, :)*phonon_field(bond_near_2(8), :)))

         !! data_13: {X_x * X_y}^2 nearest square interaction
         do i = 1 , 4
            ml_data_sub(13) = ml_data_sub(13) + sum(phonon_field(site, :)**2 * phonon_field(bond_near_1(i), :)**2)
            ml_data_sub(13) = ml_data_sub(13) + sum(phonon_field(site, :)**2 * phonon_field(bond_near_2(i), :)**2)
         end do
      end do
      ml_data_sub(9:13) = ml_data_sub(9:13)/Ns/ntime/2

      ! data_14: abs(X) average
      ml_data_sub(14) = sum(abs(phonon_field))/Ns/ntime

      ! data_15: NN time correlation
      do tau = 1, ntime
         ml_data_sub(15) = ml_data_sub(15) + sum((phonon_field(:, tau)*phonon_field(:, modi(tau + 2,ntime))))
         ml_data_sub(16) = ml_data_sub(16) + sum((phonon_field(:, tau)*phonon_field(:, ntime-tau+1)))/2d0
      end do
      
   END subroutine

   subroutine local_update_ML()
      implicit none
      integer :: time, site, bond, local_accept, sweep, tau
      real(8) :: old_field(bond_number, ntime)
      real(8) :: old_ln_cw
      real(8) :: pf_jump(1), d_energy, prob, accept_ratio
      integer :: time_list(1)
      real(8) :: distance
      old_field = phonon_field
      old_ln_cw = ln_cw
      ln_cw = 0d0
      local_accept = 0
      debug_accept = 0
      cumu_en = 0d0
      accept_ratio = 0d0
      !distance = jump_distance/abs(ep_parameter)
      distance = jump_distance
      if (.not. allocated(parameters)) call init_ml_parameter()

      ! local update
      do sweep = 1, ML_sweep
         do time = 1, ntime
            time_list(1) = time
            do bond = 1, bond_number

               pf_jump(1) = distance*drand_sym()
               call cal_phdiff_single(bond, time_list, 1, pf_jump, d_energy , spatial=.false.)
               call cal_prob_ML(prob, d_energy)
               !print*,'ML_local_prob:',prob
               if (prob > drand()) then
                  !accept
                  local_accept = local_accept + 1
                  phonon_field(bond, time_list) = phonon_field(bond, time_list) + pf_jump(1)
                  cumu_en = cumu_en + d_energy
               end if
            end do
         end do
      end do

      !print *, 'KE_dom_ratio:', real(debug_accept)/bond_number/ntime/ML_sweep
      !print *, 'accept_ratio:', real(local_accept)/bond_number/ntime/ML_sweep
      accept_ratio = real(local_accept)/bond_number/ntime/ML_sweep
      
      !! decide whether to accept this ML update or not
      d_energy = 0d0
      d_energy = d_energy + cal_phdiff_all(old_field, phonon_field, effective=.false., wolff=.false., spatial=.false.)
      call cal_new_weight()
      !print *, 'ph_diff:', d_energy
      d_energy = d_energy + (ln_cw - old_ln_cw)/(-delt)
      !print *, '(ln_cw - old_ln_cw)/(-delt):', (ln_cw - old_ln_cw)/(-delt)
      !cumu_en = cal_phdiff_all(old_field, phonon_field, effective=.true., wolff=.false., spatial=.false.)
      d_energy = d_energy - cumu_en
      !print *, 'cumu_en:', cumu_en
      !print *, 'final_d_en:', d_energy
      !print *, ''
      call cal_prob_ML(prob, d_energy)
      if (prob > drand()) then
         !accept
         ML_accept = ML_accept + 1
         ML_accept_ratio = ML_accept_ratio + accept_ratio
         !call init_g_T0()
      else
         !reject
         ML_reject = ML_reject + 1
         phonon_field = old_field
         ln_cw = old_ln_cw
      end if

   end subroutine

   !> cal the en of bond at time_indexes in time_list,time_length >= 1
   subroutine cal_phdiff_single(bond, time_list,time_length, pf_jump, energy,spatial)

      implicit none
      integer, intent(in) :: bond,time_length
      real(8) :: pf_jump(time_length), energy
      real(8) ::  E_1, E_2, ph_be, ph_af, ph_1, ph_2, ph_1_new, ph_2_new
      integer :: time_list(time_length),time
      integer ::  tau, near_bond(8)
      real(8) :: variable(n_variable)
      logical :: spatial
      energy = 0d0
      variable = 0d0
      if(time_length == 1) time = time_list(1)
      if(.not.spatial) then
         !!time direction

         ph_be = phonon_field(bond, modi(time_list(1) - 1, ntime))
         ph_af = phonon_field(bond, modi(time_list(time_length) + 1, ntime))
         ph_1 = phonon_field(bond, time_list(1))
         ph_1_new = ph_1 + pf_jump(1)
         ph_2 = phonon_field(bond, modi(time_list(time_length), ntime))
         ph_2_new = ph_2 + pf_jump(time_length)

         ! modified original ph kenetic energy
         E_1 = ((ph_be - ph_1)**2 + (ph_af - ph_2)**2)
         E_2 = ((ph_be - ph_1_new)**2 + (ph_af - ph_2_new)**2)
         variable(1) = variable(1) + (E_2 - E_1)
      end if
      !! spatial direction

      ! onsite squared term
      E_1 = sum(phonon_field(bond, time_list)**2)
      E_2 = sum((phonon_field(bond, time_list) + pf_jump)**2)
      variable(2) = variable(2) + (E_2 - E_1)

      ! onsite fourth power term
      E_1 = sum(phonon_field(bond, time_list)**4)
      E_2 = sum((phonon_field(bond, time_list) + pf_jump)**4)
      variable(3) = variable(3) + (E_2 - E_1)

      ! onsite sixth power term
      E_1 = sum(phonon_field(bond, time_list)**6)
      E_2 = sum((phonon_field(bond, time_list) + pf_jump)**6)
      variable(4) = variable(4) + (E_2 - E_1)

      ! nearest and next nearest interaction
      call get_near_bond(bond, near_bond)

      ! nearest anti interaction
      
      E_2 = sum((matmul(pf_jump,transpose(phonon_field(near_bond(1:2), time_list)))))
      variable(5) = variable(5) + (E_2)

      !  nearest ferro interaction
      
      E_2 = sum((matmul(pf_jump,transpose(phonon_field(near_bond(3:4), time_list)))))
      variable(6) = variable(6) + (E_2)

      !  next nearest longitude interaction
      
      E_2 = sum((matmul(pf_jump,transpose(phonon_field(near_bond(5:6), time_list)))))
      variable(7) = variable(7) + (E_2)

      ! next nearest transverse interaction
      
      E_2 = sum((matmul(pf_jump,transpose(phonon_field(near_bond(7:8), time_list)))))
      variable(8) = variable(8) + (E_2)

      ! nearest {X_{x} * X_{y}}^2
      E_1 = sum(matmul(phonon_field(bond, time_list)**2,transpose(phonon_field(near_bond(1:4), time_list)**2)))
      E_2 = sum(matmul((phonon_field(bond, time_list) + pf_jump)**2, transpose(phonon_field(near_bond(1:4), time_list)**2)))
      variable(9) = variable(9) + (E_2 - E_1)

      ! abs(X) average
      E_1 = sum(abs(phonon_field(bond, time_list)))
      E_2 = sum(abs(phonon_field(bond, time_list) + pf_jump))
      variable(10) = variable(10) + (E_2 - E_1)

      energy = sum(variable*parameters)

      ! ph kenetic energy
      energy = energy + variable(1)*(0.5/(delt*omega)**2)
      ! ph potential energy
      energy = energy + variable(2) * (0.5d0*D)

   end subroutine

   subroutine cal_prob_ML(p, energy)
      implicit none
      real(8) :: p, energy
      real(8) :: R
      R = exp(-delt*energy)
      if (R >= 0d0) then
         if (Metropolis) then
            p = real(R)
            if (greatest_decent) then
               p = floor(p)
            end if
         else!heat bath
            p = real(R/(1 + R)) ! accept probobility
         end if
         !print*,"P = ",prob

      else
         print *, "SIGN PROBLEM,R:", R, 'in ML'
         stop
      end if

   end subroutine

   real(8) function cal_phdiff_all(old_field, new_field, effective, wolff, spatial) result(en)
      implicit none
      real(8), intent(in) :: old_field(bond_number, ntime), new_field(bond_number, ntime)
      real(8) :: E_1, E_2
      integer :: tau, site, bond ,bond_near_1(8) , bond_near_2(8), i
      logical :: effective, wolff, spatial
      real(8) :: variable(n_variable)
      real(8) :: old_field_bond(ntime), new_field_bond(ntime)
      if (.not. allocated(parameters)) call init_ml_parameter()
      en = 0d0
      variable = 0d0

      if (.not. wolff) then


         !print*,'potential en:',en
         if (.not. spatial) then
            do tau = 1, ntime

               E_1 = sum((old_field(:, tau) - old_field(:, modi(tau + 1, ntime)))**2)
               E_2 = sum((new_field(:, tau) - new_field(:, modi(tau + 1, ntime)))**2)
               variable(1) = variable(1) + E_2 - E_1
               !    print*,'kenetic en:', E_2 - E_1
            end do
         end if

         variable(2) = sum(new_field**2 - old_field**2)
      end if

      if (effective) then ! effective energy difference
         if(.not.wolff) then
         variable(3) = sum(new_field**4 - old_field**4)
         variable(4) = sum(new_field**6 - old_field**6)
         end if
         E_1 = 0d0
         E_2 = 0d0
         do site = 1, Ns
            call get_near_bond(site, bond_near_1)
            call get_near_bond(site + Ns, bond_near_2)

            !!var5 : nearest anti interaction

            E_1 = E_1 + sum(old_field(site, :)*old_field(bond_near_1(1), :))
            E_1 = E_1 + sum(old_field(site, :)*old_field(bond_near_1(2), :))
            E_1 = E_1 + sum(old_field(site+Ns, :)*old_field(bond_near_2(1), :))
            E_1 = E_1 + sum(old_field(site+Ns, :)*old_field(bond_near_2(2), :))

            E_2 = E_2 + sum(new_field(site, :)*new_field(bond_near_1(1), :))
            E_2 = E_2 + sum(new_field(site, :)*new_field(bond_near_1(2), :))
            E_2 = E_2 + sum(new_field(site+Ns, :)*new_field(bond_near_2(1), :))
            E_2 = E_2 + sum(new_field(site+Ns, :)*new_field(bond_near_2(2), :))

            variable(5) = variable(5) + (E_2 - E_1)
            E_1 = 0d0
            E_2 = 0d0
            !!var6 : nearest ferro interaction

            E_1 = E_1 + sum(old_field(site, :)*old_field(bond_near_1(3), :))
            E_1 = E_1 + sum(old_field(site, :)*old_field(bond_near_1(4), :))
            E_1 = E_1 + sum(old_field(site+Ns, :)*old_field(bond_near_2(3), :))
            E_1 = E_1 + sum(old_field(site+Ns, :)*old_field(bond_near_2(4), :))

            E_2 = E_2 + sum(new_field(site, :)*new_field(bond_near_1(3), :))
            E_2 = E_2 + sum(new_field(site, :)*new_field(bond_near_1(4), :))
            E_2 = E_2 + sum(new_field(site+Ns, :)*new_field(bond_near_2(3), :))
            E_2 = E_2 + sum(new_field(site+Ns, :)*new_field(bond_near_2(4), :))

            variable(6) = variable(6) + (E_2 - E_1)
            E_1 = 0d0
            E_2 = 0d0
            !!var7 : X_{i,l} * X_{i+1,l} next nearest long

            !La direction
            E_1 = E_1 + sum((old_field(site, :)*old_field(bond_near_1(5), :)))
            E_1 = E_1 + sum((old_field(site, :)*old_field(bond_near_1(6), :)))

            E_2 = E_2 + sum((new_field(site, :)*new_field(bond_near_1(5), :)))
            E_2 = E_2 + sum((new_field(site, :)*new_field(bond_near_1(6), :)))
            !Lb direction
            E_1 = E_1 + sum((old_field(site+Ns, :)*old_field(bond_near_2(5), :)))
            E_1 = E_1 + sum((old_field(site+Ns, :)*old_field(bond_near_2(6), :)))
            E_2 = E_2 + sum((new_field(site+Ns, :)*new_field(bond_near_2(5), :)))
            E_2 = E_2 + sum((new_field(site+Ns, :)*new_field(bond_near_2(6), :)))

            variable(7) = variable(7) + (E_2 - E_1)
            E_1 = 0d0
            E_2 = 0d0
            !!var8 :X_{i,l} * X_{i+1,l} next nearest trans interaction

            !La direction
            E_1 = E_1 + sum((old_field(site, :)*old_field(bond_near_1(7), :)))
            E_1 = E_1 + sum((old_field(site, :)*old_field(bond_near_1(8), :)))
            E_2 = E_2 + sum((new_field(site, :)*new_field(bond_near_1(7), :)))
            E_2 = E_2 + sum((new_field(site, :)*new_field(bond_near_1(8), :)))
            !Lb direction
            E_1 = E_1 + sum((old_field(site+Ns, :)*old_field(bond_near_2(7), :)))
            E_1 = E_1 + sum((old_field(site+Ns, :)*old_field(bond_near_2(8), :)))
            E_2 = E_2 + sum((new_field(site+Ns, :)*new_field(bond_near_2(7), :)))
            E_2 = E_2 + sum((new_field(site+Ns, :)*new_field(bond_near_2(8), :)))
            variable(8) = variable(8) + (E_2 - E_1)
            E_1 = 0d0
            E_2 = 0d0
            !! var9 : {X_x * X_y}^2 nearest square interaction
            do i = 1 , 4
               E_1 = E_1 + sum(old_field(site, :)**2 * old_field(bond_near_1(i), :)**2)
               E_1 = E_1 + sum(old_field(site+Ns, :)**2 * old_field(bond_near_2(i), :)**2)
               E_2 = E_2 + sum(new_field(site, :)**2 * new_field(bond_near_1(i), :)**2)
               E_2 = E_2 + sum(new_field(site+Ns, :)**2 * new_field(bond_near_2(i), :)**2)
            end do
            variable(9) = variable(9) + (E_2 - E_1)
            E_1 = 0d0
            E_2 = 0d0
         end do
         variable(5:9) = variable(5:9)/2

         variable(10) = sum(abs(new_field) - abs(old_field))


      end if
      en = variable(1) * (0.5/(delt*omega)**2) + variable(2) * (0.5d0*D)
      if(effective) en = en + sum(variable*parameters)

   end function
   subroutine spatial_domain_update_ML()
      implicit none
      integer :: Kx, Ky,K_list(2,2),local_accept,local_count,local_accept_2pi
      real(8) :: phi_x, phi_y
      real(8) :: pf_jump(bond_number) , pf_temp(bond_number, ntime)
      real(8) :: d_energy
      integer :: bond, time, tau, bond_1, i
      real(8) :: old_field(bond_number, ntime), old_ln_cw
      integer :: spatial_sweep = 1, sweep_count
      real(8) :: accept_ratio
      old_field = phonon_field
      old_ln_cw = ln_cw
      ln_cw = 0d0
      cumu_en = 0d0
      local_accept = 0
      local_count = 0
      local_accept_2pi = 0
      accept_ratio = 0d0
      !K_list = reshape((/ La/2+1, 1, La/2+2, Lb/2+1, 1/), (/2, 2/))
      do sweep_count = 1 , spatial_sweep
         do Kx = 1 , La
            do Ky = 1 , Lb
               !Kx = irand(La) + 1
               !Ky = irand(Lb) + 1
               call generate_newfield_space(pf_jump, irand(La)+1, irand(Lb)+1,& 
               & ml_distance_ratio*global_update_distance/abs(ep_parameter))
               !print*,'pf_jump(k):',pf_jump
               do bond = 1 , bond_number
                  pf_temp(bond, :) = phonon_field(bond,:) + pf_jump(bond)
               end do

               d_energy = cal_phdiff_all(phonon_field, pf_temp, effective=.true., wolff=.false., spatial=.true.)
               call cal_prob_ML(prob, d_energy)

               if (prob > drand()) then
                  !accept
                  phonon_field = pf_temp
                  cumu_en = cumu_en + d_energy
                  local_accept = local_accept + 1
                  !print*, 'accept: kx,ky:', Kx, Ky
               end if

               call generate_newfield_space(pf_jump, La/2+1, Lb/2+1, & 
               & ml_distance_ratio*global_update_distance/abs(ep_parameter))
               !print*,'pf_jump(k):',pf_jump
               do bond = 1 , bond_number
                  pf_temp(bond, :) = phonon_field(bond,:) + pf_jump(bond)
               end do
               d_energy = cal_phdiff_all(phonon_field, pf_temp, effective=.true., wolff=.false., spatial=.true.)
               call cal_prob_ML(prob, d_energy)

               if (prob > drand()) then
                  !accept
                  phonon_field = pf_temp
                  cumu_en = cumu_en + d_energy
                  local_accept_2pi = local_accept_2pi + 1
                  !print*, 'accept: pi,pi'
               end if
               local_count = local_count + 1
            end do
         end do
      end do
      ! print *, 'accept_ratio:', real(local_accept)/local_count
      !print *, 'accept_ratio_2pi:', real(local_accept_2pi)/local_count
      accept_ratio = real(local_accept + local_accept_2pi)/local_count / 2
      ! decide whether to accept this kspace update or not
      d_energy = 0d0
      d_energy = d_energy + cal_phdiff_all(old_field, phonon_field, effective=.false., wolff=.false., spatial=.true.)
      call cal_new_weight()
      !print *, 'ph_diff:', d_energy
      d_energy = d_energy + (ln_cw - old_ln_cw)/(-delt)
      !print *, '(ln_cw - old_ln_cw)/(-delt):', (ln_cw - old_ln_cw)/(-delt)
      !cumu_en = cal_phdiff_all(old_field, phonon_field, effective=.true., wolff=.false., spatial=.true.)
      d_energy = d_energy - cumu_en
      !print *, 'cumu_en:', cumu_en
      !print *, 'final_d_en:', d_energy
      !print *, ''
      call cal_prob_ML(prob, d_energy)
      if (prob > drand()) then
         !accept
         ML_s_accept = ML_s_accept + 1
         !call init_g_T0()
         ml_s_accept_ratio = ml_s_accept_ratio + accept_ratio
      else
         !reject
         ML_s_reject = ML_s_reject + 1
         phonon_field = old_field
         ln_cw = old_ln_cw
      end if

   end subroutine
   subroutine time_domain_update_ML()
      implicit none
      integer :: k_index,k_temp
      real(8) :: phi
      real(8) :: pf_jump(ntime)
      real(8) :: d_energy, temp_en
      integer :: bond, time, tau, bond_1, bond_2,local_accept
      real(8) :: old_field(bond_number, ntime), old_ln_cw
      REAL(8) :: accept_ratio 
      old_field = phonon_field
      old_ln_cw = ln_cw
      ln_cw = 0d0
      cumu_en = 0d0
      d_energy = 0d0
      local_accept = 0
      n_k_ML = 20
      if(.not. allocated(parameters)) call init_ml_parameter()
      do bond = 1, bond_number
         do k_index = 1 , n_k_ML
            k_temp = irand(ntime/2) + 1
            !k_temp = k_index
            call generate_newfield_time(pf_jump, k_temp, jump_distance)
            !print*,'pf_jump(k):',pf_jump
            !do tau = 1, ntime
            !    call cal_phdiff_single(bond, tau, pf_jump(tau), temp_en)
            !    phonon_field(bond, tau) = phonon_field(bond, tau) + pf_jump(tau)
            !    d_energy = d_energy + temp_en
            !end do
            call cal_phdiff_time(bond, pf_jump, d_energy)
            call cal_prob_ML(prob, d_energy)

            if (prob > drand()) then
               !accept
               !print*, 'accept: k:', k_temp
               phonon_field(bond, :) = phonon_field(bond, :) + pf_jump
               cumu_en = cumu_en + d_energy
               local_accept = local_accept + 1
            else
               !print*, 'reject: k:', k_temp
               ! phonon_field(bond, :) = phonon_field(bond, :) - pf_jump
            end if
            d_energy = 0d0
         end do

      end do
      
      !print*, 'accept_ratio:', real(local_accept)/bond_number/n_k_ML
      accept_ratio = real(local_accept)/bond_number/n_k_ML
      ! decide whether to accept this kspace update or not
      d_energy = 0d0
      d_energy = d_energy + cal_phdiff_all(old_field, phonon_field, effective=.false., wolff=.false., spatial=.false.)
      call cal_new_weight()
      !print *, 'ph_diff:', d_energy
      d_energy = d_energy + (ln_cw - old_ln_cw)/(-delt)
      !print *, '(ln_cw - old_ln_cw)/(-delt):', (ln_cw - old_ln_cw)/(-delt)
      !cumu_en = cal_phdiff_all(old_field, phonon_field, effective=.true., wolff=.false., spatial=.false.)
      d_energy = d_energy - cumu_en
      !print *, 'cumu_en:', cumu_en
      !print *, 'final_d_en:', d_energy
      !print *, ''
      call cal_prob_ML(prob, d_energy)
      if (prob > drand()) then
         !accept
         ML_t_accept = ML_t_accept + 1
         ML_t_accept_ratio = ML_t_accept_ratio + accept_ratio
         !call init_g_T0()
      else
         !reject
         ML_t_reject = ML_t_reject + 1
         phonon_field = old_field
         ln_cw = old_ln_cw
      end if

   end subroutine

   subroutine cal_phdiff_time(bond, pf_jump, en)
      implicit none
      integer, intent(in) :: bond
      real(8), intent(in) :: pf_jump(ntime)
      real(8) :: en
      real(8) :: E_1, E_2, variable(n_variable)
      real(8) :: new_field(ntime), old_field(ntime)
      integer :: tau, near_bond(8)
      en = 0d0
      variable = 0d0

      new_field = phonon_field(bond, :) + pf_jump
      old_field = phonon_field(bond, :)

      !print*,'potential en:',en
      do tau = 1, ntime

         E_1 = (old_field(tau) - old_field(modi(tau + 1, ntime)))**2
         E_2 = (new_field(tau) - new_field(modi(tau + 1, ntime)))**2

         variable(1) = variable(1) + E_2 - E_1
         !    print*,'kenetic en:', E_2 - E_1
      end do
      !print*, 'K_en:', variable(1) * (0.5/(delt*omega)**2)
      !stop
      ! onsite squared term
      variable(2) = sum(new_field**2 - old_field**2)
      variable(3) = sum(new_field**4 - old_field**4)
      variable(4) = sum(new_field**6 - old_field**6)
      call get_near_bond(bond, near_bond)
      ! nearest anti interaction
      E_1 = sum(matmul(phonon_field(near_bond(1:2), :),old_field))
      E_2 = sum(matmul(phonon_field(near_bond(1:2), :),new_field))
      variable(5) = variable(5) + (E_2 - E_1)

      !  nearest ferro interaction
      E_1 = sum(matmul(phonon_field(near_bond(3:4), :),old_field))
      E_2 = sum(matmul(phonon_field(near_bond(3:4), :),new_field))
      variable(6) = variable(6) + (E_2 - E_1)

      !  next nearest longitude interaction
      E_1 = sum(matmul(phonon_field(near_bond(5:6), :),old_field))
      E_2 = sum(matmul(phonon_field(near_bond(5:6), :),new_field))
      variable(7) = variable(7) + (E_2 - E_1)

      ! next nearest transverse interaction
      E_1 = sum(matmul(phonon_field(near_bond(7:8), :),old_field))
      E_2 = sum(matmul(phonon_field(near_bond(7:8), :),new_field))
      variable(8) = variable(8) + (E_2 - E_1)

      ! nearest {X_{x} * X_{y}}^2
      E_1 = sum(matmul(phonon_field(near_bond(1:4), :)**2,old_field**2))
      E_2 = sum(matmul(phonon_field(near_bond(1:4), :)**2,(new_field)**2))
      variable(9) = variable(9) + (E_2 - E_1)

      ! abs(X) average
      E_1 = sum(abs(old_field))
      E_2 = sum(abs(new_field))
      variable(10) = variable(10) + (E_2 - E_1)



      en = sum(variable*parameters)
      en = en + variable(1)*(0.5/(delt*omega)**2) + variable(2)*0.5d0*D



   end subroutine
   subroutine Wolff_update()

      implicit none
      integer :: bond, start_bond

      real(8) :: prob, d_energy
      start_bond = irand(bond_number) + 1

      call init_wait_list(bond_number)
      call init_cluster(bond_number)

      !start from a random bond
      bond = start_bond
      d_energy = 1000d0 ! large enough to accept the first bond
      call add_wait_list(bond, d_energy)

      do while (current_wl > 0)
         ! pop a bond from wait_list and decide whether to add it into cluster
         bond = wait_list(current_wl)%bond
         d_energy = wait_list(current_wl)%en
         current_wl = current_wl - 1
         !print*, 'd_energy:', d_energy
         call cal_prob_ML(prob, d_energy)
         prob = 1d0 - prob
         !print*,'prob:',prob
         if (prob > drand()) then
            !print*,'accept:',bond
            call add_cluster(bond)
            call add_nearest_bond(bond)
         end if
      end do
      !print*,'current_cluster_x:',current_cluster
      if (current_cluster /= bond_number) then
         call cluster_flip()
      else
         !call cluster_flip()
          wolff_reject = wolff_reject + 1
      end if
      deallocate(cluster, wait_list)

   end subroutine
   subroutine time_wolff_update()
      implicit none
      integer :: tau,start_tau
      real(8) :: d_energy
      real(8) :: old_field(bond_number, ntime), old_ln_cw
      old_field = phonon_field
      old_ln_cw = ln_cw
      ln_cw = 0d0
      cumu_en = 0d0
      start_tau = irand(ntime) + 1

      call init_wait_list(2)
      call init_cluster(ntime)

      !start from a random time

      tau = start_tau
      d_energy = 1000d0 ! large enough to accept the first time
      call add_wait_list(tau, d_energy)

      do while (current_wl > 0)
         ! pop a time slice from wait_list and decide whether to add it into cluster
         tau = wait_list(current_wl)%bond
         d_energy = wait_list(current_wl)%en
         current_wl = current_wl - 1
         !print*, 'd_energy:', d_energy
         call cal_prob_ML(prob, d_energy)
         prob = 1d0 - prob
         !print*,'prob:',prob
         if (prob > drand()) then
            !print*,'accept:',bond
            call add_cluster(tau)
            call add_nearest_time(tau)
         else
            !cumu_en = cumu_en + d_energy
         end if
      end do
      if(current_cluster /= ntime) then
         call cluster_time_flip()
      else
         wolff_time_reject = wolff_time_reject + 1
      end if
      deallocate(cluster, wait_list)
   end subroutine
   subroutine add_nearest_time(time)
      implicit none
      integer :: time
      integer :: tau
      real(8) :: d_energy

      tau = modi(time - 1, ntime)
      if(findloc(cluster,tau,1) == 0) then
         d_energy = sum(4 * phonon_field(:, tau) * phonon_field(:, time) * (0.5/(delt*omega)**2))
         call add_wait_list(tau, d_energy)
      end if

      tau = modi(time + 1, ntime)
      if(findloc(cluster,tau,1) == 0) then
         d_energy = sum(4 * phonon_field(:, tau) * phonon_field(:, time) * (0.5/(delt*omega)**2))
         call add_wait_list(tau, d_energy)
      end if
   end subroutine
   subroutine init_wait_list(n)
      implicit none
      integer :: n
      integer :: i
      if (.not. allocated(wait_list)) allocate (wait_list(n))
      do i = 1, n
         wait_list(i)%bond = 0
         wait_list(i)%en = 0
      end do
      current_wl = 0
   end subroutine

   subroutine init_cluster(n)
      implicit none
      integer :: n
      integer :: i
      if (.not. allocated(cluster)) allocate (cluster(n))
      do i = 1, n
         cluster(i) = 0
      end do
      current_cluster = 0
   end subroutine

   subroutine add_wait_list(bond, d_energy)
      implicit none
      integer :: bond, bond_ind
      real(8) :: d_energy


      bond_ind = findloc(wait_list(1:current_wl)%bond, bond, 1)
      if (bond_ind /= 0) then
         wait_list(bond_ind)%en = wait_list(bond_ind)%en + d_energy
         return
      else
         if (current_wl == size(wait_list)) stop 'wait_list is full'
         current_wl = current_wl + 1
         wait_list(current_wl)%bond = bond
         wait_list(current_wl)%en = d_energy
      end if
   end subroutine

   subroutine add_cluster(bond)
      implicit none
      integer :: bond
      if (current_cluster == size(cluster)) then
         print *, 'cluster is full'
         print *, 'bond:', bond
         print *, 'cluster:', cluster
      end if
      current_cluster = current_cluster + 1
      cluster(current_cluster) = bond

   end subroutine

   subroutine add_nearest_bond(bond)
      !> add the nearest bonds into wait_list
      implicit none
      integer :: bond
      integer :: bond_nearest(8)
      integer :: i
      real(8) :: d_energy

      call get_near_bond(bond, bond_nearest)
      !print*,'bond:',bond
      do i = 1, 8
         if (bond_nearest(i) > bond_number) then
            print *, 'bond_nearest(i) > bond_number', bond_nearest
            print *, 'bond', bond
         end if
         !print*, 'bond_nearest(i):', bond_nearest(i)
         if (findloc(cluster(1:current_cluster), bond_nearest(i), 1) == 0) then
            if (i <= 2) then ! transverse field
               d_energy = -2*J_anti*sum(phonon_field(bond, :)*phonon_field(bond_nearest(i), :))
            elseif (i <= 4) then ! longitude field
               d_energy = -2*J_ferro*sum(phonon_field(bond, :)*phonon_field(bond_nearest(i), :))
            elseif (i <= 6) then
               d_energy = -2*J_long*sum(phonon_field(bond, :)*phonon_field(bond_nearest(i), :))
            elseif (i <= 8) then
               d_energy = -2*J_trans*sum(phonon_field(bond, :)*phonon_field(bond_nearest(i), :))
            end if

            !print*,'d_energy:',d_energy
         call add_wait_list(bond_nearest(i), d_energy)
         end if
         

      end do
      
   end subroutine

   subroutine get_near_bond(bond, bond_nearest)
      !> get the nearest bond of bond
      !> 1,2 : nearest_anti , 3,4 : nearest_ferro, 5,6 : next_nearest_long, 7,8 : next_nearest_trans
      implicit none
      integer,intent(in) :: bond
      integer :: site,site_up,site_down,site_left,site_right,site_downright,site_upleft
      integer :: bond_nearest(8)
      integer :: direction,bond_1
      bond_nearest = 0
      direction = bond_list(bond)%direction
      site = bond_list(bond)%site_1
      site_up = modi(site+La, Ns)
      site_down = modi(site-La, Ns)
      site_left = site - 1
      if (modi(site_left, La) == La) site_left = site_left + La
      site_right = site + 1
      if (modi(site_right, La) == 1) site_right = site_right - La
      site_downright = site_down + 1
      if (modi(site_downright, La) == 1) site_downright = site_downright - La
      site_upleft = site_up - 1
      if (modi(site_upleft, La) == La) site_upleft = site_upleft + La

      if (direction == 1) then
         ! x direction bond
         ! nearest_anti
         bond_nearest(1) = bond + Ns
         bond_nearest(2) = site_downright + Ns
         ! nearest_ferro
         bond_nearest(3) = site_down + Ns
         bond_nearest(4) = site_right + Ns
         ! next_nearest_long
         bond_1 = bond + 1
         if (modi(bond_1, La) == 1) bond_1 = bond_1 - La
         bond_nearest(5) = bond_1
         bond_1 = bond - 1
         if (modi(bond_1, La) == La) bond_1 = bond_1 + La
         bond_nearest(6) = bond_1
         ! next_nearest_trans
         bond_1 = bond + La
         if (bond_1 > Ns) bond_1 = bond_1 - Ns
         bond_nearest(7) = bond_1
         bond_1 = bond - La
         if (bond_1 <= 0) bond_1 = bond_1 + Ns
         bond_nearest(8) = bond_1

      elseif (direction == 2) then
         ! y direction bond
         ! nearest_anti
         bond_nearest(1) = bond - Ns
         bond_nearest(2) = site_upleft
         ! nearest_ferro
         bond_nearest(3) = site_up
         bond_nearest(4) = site_left

         ! next_nearest_long
         bond_1 = bond + La
         if (bond_1 > 2*Ns) bond_1 = bond_1 - Ns
         bond_nearest(5) = bond_1
         bond_1 = bond - La
         if (bond_1 <= Ns) bond_1 = bond_1 + Ns
         bond_nearest(6) = bond_1
         ! next_nearest_trans
         bond_1 = bond + 1
         if (modi(bond_1, La) == 1) bond_1 = bond_1 - La
         bond_nearest(7) = bond_1
         bond_1 = bond - 1
         if (modi(bond_1, La) == La) bond_1 = bond_1 + La
         bond_nearest(8) = bond_1
      end if

   end subroutine

   subroutine cluster_flip()

      implicit none

      real(8) :: old_field(bond_number, ntime)
      real(8) :: old_ln_cw
      real(8) :: d_energy
      old_field = phonon_field
      old_ln_cw = ln_cw
      d_energy = 0d0
      ln_cw = 0d0
      ! calculate the total energy difference
     ! print *, 'current_cluster:', current_cluster
      !print*,'cluster:',cluster(1:current_cluster)
      phonon_field(cluster(1:current_cluster), :) = -phonon_field(cluster(1:current_cluster), :)
      !d_energy = d_energy + cal_phdiff_all(old_field, phonon_field, effective=.false.)
      call cal_new_weight()
      !print *, 'ph_diff:', d_energy
      d_energy = d_energy + (ln_cw - old_ln_cw)/(-delt)
      !print *, '(ln_cw - old_ln_cw)/(-delt):',d_energy
      cumu_en = cal_phdiff_all(old_field, phonon_field, effective=.true., wolff=.true., spatial=.true.)
      d_energy = d_energy - cumu_en
      !print *, 'cumu_en:', cumu_en
      !print *, 'final_d_en:', d_energy
      
      call cal_prob_ML(prob, d_energy)
      if (prob > drand()) then
         !accept
         Wolff_accept = wolff_accept + 1
         !print *, 'current_cluster:', current_cluster
         !print *, ''
         !call init_g_T0()
      else
         !reject
         WolFF_reject = wolff_reject + 1
         phonon_field = old_field
         ln_cw = old_ln_cw
      end if
   end subroutine

   subroutine cluster_time_flip()

      implicit none

      real(8) :: old_field(bond_number, ntime)
      real(8) :: old_ln_cw
      real(8) :: d_energy
      old_field = phonon_field
      old_ln_cw = ln_cw
      d_energy = 0d0
      ln_cw = 0d0
      ! calculate the total energy difference
      print *, 'current_cluster:', current_cluster
      !print*,'cluster:',cluster(1:current_cluster)
      phonon_field(:, cluster(1:current_cluster)) = - phonon_field(:, cluster(1:current_cluster))
      !d_energy = d_energy + cal_phdiff_all(old_field, phonon_field, effective=.false.)
      call cal_new_weight()
      !print *, 'ph_diff:', d_energy
      d_energy = d_energy + (ln_cw - old_ln_cw)/(-delt)
      !print *, '(ln_cw - old_ln_cw)/(-delt):', (ln_cw - old_ln_cw)/(-delt)
      !cumu_en = cal_phdiff_all(old_field, phonon_field, effective=.true., wolff=.true., spatial=.true.)
      !d_energy = d_energy - cumu_en
      !print *, 'cumu_en:', cumu_en
      print *, 'final_d_en:', d_energy
      !print *, ''
      call cal_prob_ML(prob, d_energy)
      if (prob > drand()) then
         !accept
         wolff_time_accept = wolff_time_accept + 1
         !call init_g_T0()
      else
         !reject
         wolff_time_reject = wolff_time_reject + 1
         phonon_field = old_field
         ln_cw = old_ln_cw
      end if
   end subroutine

end module
