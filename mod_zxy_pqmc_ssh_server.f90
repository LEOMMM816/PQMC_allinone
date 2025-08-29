MODULE ZXYP_SSH_server
  !use mpi
  use pf_setting
  use ran_state
  IMPLICIT NONE

#ifdef MPI
  include 'mpif.h'
#endif

CONTAINS


  
 
  complex(8) function get_bond_coupling(bond,time,phase) result(re)

    implicit none
    integer,intent(in) :: bond,time
    logical :: phase
    if(phase) then
      !re = (-hop) * bond_list(bond)%phase * exp(complex(0d0,TBC_phase_number((bond-1)/Ns + 1)*pi /La))
    else
      !re = (-hop)* exp(complex(0d0,TBC_phase_number((bond-1)/Ns + 1)*pi /La))
    end if
  end function

  subroutine output_phonon_field()
    implicit none
    integer :: time
    character*4 :: ci3
    character*40 :: str = 'boson_field.dat'
#ifdef MPI
    write (ci3, '(1i4)') myid + 1
    str = 'ph_field_core'//trim(adjustl(ci3))//'.dat'
    !if (mod(myid, MPI_one_block) /= 0) return
#endif
    write (ci3, '(1i4)') Ns
    !inquire(file = str,exist = file_exist)
    !if (.not.file_exist) then
    !   open(30,file = str,status = 'new')
    !   close(30)
    !end if
104 open (30 + myid, file=str, status='replace', err=103)
    do time = 1, ntime
      write (30 + myid, '('//trim(adjustl(ci3))//'f18.6)', err=104) boson_field(:, time)
    end do
    close (30)
103 time = time + 1
  end subroutine

!> doesn't recalculate ln_cw

end MODULE
