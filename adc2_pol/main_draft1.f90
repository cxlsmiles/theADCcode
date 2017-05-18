program main
  
  use constants
  use parameters
  use read_param
  
  implicit none
  
  
  integer :: i,j
  real(d) :: time

  integer, dimension(2) :: shp

!!$ Setting up the multiplication table

  shp(:)=(/ 8,8 /)
  MT=reshape(mtrow, shp)

!!$ Reading basic information from Molcas generated files

  call read_molcas0()
  
  allocate(e(nBas),occNum(nBas),orbSym(nBas),roccnum(nBas))
  allocate(x_dipole(nBas,nBas),y_dipole(nBas,nBas),z_dipole(nBas,nBas),dpl(nBas,nBas))
  
  call read_user()
  call read_molcas1(info)
  

!!$ Reading user's data

!  call read_user()
  call check_user()

!!$ Rearranging orbitals such that occ. orbs preceed unoccupied orbs.

  call rearrange_occ()

  if (tranmom .eq. 'x') then
     dpl(:,:)=x_dipole(:,:)
     tranflag='y'
  elseif (tranmom .eq. 'y') then
     dpl(:,:)=y_dipole(:,:)
      tranflag='y'
  elseif (tranmom .eq. 'z') then
     dpl(:,:)=z_dipole(:,:)
     tranflag='y'
  end if
     

  if (method .eq. 0) then
     write(6,*) "Activating Fano"
     call master_fano()
  elseif (method .eq. 1) then
     write(6,*) "Activating TDA"
     call master_tda()
  elseif (method .eq. 2) then
     write(6,*) "Activating ADC2"
     call master_adc2()
  elseif (method .eq. 3) then
     write(6,*) "Activating ADC2 extended"
     call master_adc2ext()
  end if

  call cpu_time(time)
  write(6,*) 'Final Time=',time," s"
  
end program main
  
  
  
