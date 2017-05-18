
module def_config
implicit none

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type Molecular_Orbital

! indice of the MO where is the hole or the particle
  integer :: ind
! spin can be a (=alpha) or b(=beta)
  character(1) :: spin

end type Molecular_Orbital

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type config

 integer :: nholes, nparts
 type(Molecular_Orbital), dimension(10) :: ihole, ipart

end type config

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type spin_adapted_config

 integer :: spin
 integer :: typ
 integer :: nconfig
 double precision, dimension(:), allocatable :: coeff
 type(config), dimension(:), allocatable :: cfg

 integer :: nholes, nparts
 integer, dimension(10) :: ihole, ipart

end type spin_adapted_config

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine def_spac(spac)

type(spin_adapted_config), intent(inout) :: spac

! For doubly-ionized states,
! spin adapted configurations are defined as Francesco Tarantelli's Paper (Chemical Physics 329 (2006) 11-21)

if(spac%spin==1) then


  if(spac%nholes==2) then
      if(spac%ihole(1) == spac%ihole(2)) then

        spac%nconfig = 1
        if(.not. allocated(spac%coeff)) allocate(spac%coeff(spac%nconfig),spac%cfg(spac%nconfig))
        spac%cfg(:)%nholes = 3
        spac%cfg(:)%nparts = 1

        spac%coeff(1) = 1d0
        spac%cfg(1)%ihole(1)%ind  = spac%ihole(1)
        spac%cfg(1)%ihole(1)%spin = 'a'
        spac%cfg(1)%ihole(2)%ind  = spac%ihole(2)
        spac%cfg(1)%ihole(2)%spin = 'b'

        spac%cfg(1)%ihole(3)%ind  = 0
        spac%cfg(1)%ihole(3)%spin = 'n'

      else

        spac%nconfig = 2
        if(.not. allocated(spac%coeff)) allocate(spac%coeff(spac%nconfig),spac%cfg(spac%nconfig))
        spac%cfg(:)%nholes = 3
        spac%cfg(:)%nparts = 1

        spac%coeff(1) = -1d0/dsqrt(2d0)
        spac%cfg(1)%ihole(1)%ind = spac%ihole(1)
        spac%cfg(1)%ihole(1)%spin = 'a'
        spac%cfg(1)%ihole(2)%ind = spac%ihole(2)
        spac%cfg(1)%ihole(2)%spin = 'b'

        spac%cfg(1)%ihole(3)%ind  = 0
        spac%cfg(1)%ihole(3)%spin = 'n'

        spac%coeff(2) = 1d0/dsqrt(2d0)
        spac%cfg(2)%ihole(1)%ind = spac%ihole(1)
        spac%cfg(2)%ihole(1)%spin = 'b'
        spac%cfg(2)%ihole(2)%ind = spac%ihole(2)
        spac%cfg(2)%ihole(2)%spin = 'a'

        spac%cfg(2)%ihole(3)%ind  = 0
        spac%cfg(2)%ihole(3)%spin = 'n'

      endif

  elseif(spac%nholes==3) then

      if(spac%typ==1) then

        spac%nconfig = 2
        if(.not. allocated(spac%coeff)) allocate(spac%coeff(spac%nconfig),spac%cfg(spac%nconfig))
        spac%cfg(:)%nholes = 3
        spac%cfg(:)%nparts = 1

        spac%coeff(1) = 1d0/dsqrt(2d0)
        spac%cfg(1)%ipart(1)%ind = spac%ipart(1)
        spac%cfg(1)%ipart(1)%spin = 'a'
        spac%cfg(1)%ihole(1)%ind = spac%ihole(1)
        spac%cfg(1)%ihole(1)%spin = 'a'
        spac%cfg(1)%ihole(2)%ind = spac%ihole(2)
        spac%cfg(1)%ihole(2)%spin = 'a'
        spac%cfg(1)%ihole(3)%ind = spac%ihole(3)
        spac%cfg(1)%ihole(3)%spin = 'b'

        spac%coeff(2) = 1d0/dsqrt(2d0)
        spac%cfg(2)%ipart(1)%ind = spac%ipart(1)
        spac%cfg(2)%ipart(1)%spin = 'b'
        spac%cfg(2)%ihole(1)%ind = spac%ihole(1)
        spac%cfg(2)%ihole(1)%spin = 'b'
        spac%cfg(2)%ihole(2)%ind = spac%ihole(2)
        spac%cfg(2)%ihole(2)%spin = 'a'
        spac%cfg(2)%ihole(3)%ind = spac%ihole(3)
        spac%cfg(2)%ihole(3)%spin = 'b'

      elseif(spac%typ==21) then

        spac%nconfig = 4
        if(.not. allocated(spac%coeff)) allocate(spac%coeff(spac%nconfig),spac%cfg(spac%nconfig))
        spac%cfg(:)%nholes = 3
        spac%cfg(:)%nparts = 1

        spac%coeff(1) = 0.5d0
        spac%cfg(1)%ipart(1)%ind = spac%ipart(1)
        spac%cfg(1)%ipart(1)%spin = 'a'
        spac%cfg(1)%ihole(1)%ind = spac%ihole(1)
        spac%cfg(1)%ihole(1)%spin = 'b'
        spac%cfg(1)%ihole(2)%ind = spac%ihole(2)
        spac%cfg(1)%ihole(2)%spin = 'a'
        spac%cfg(1)%ihole(3)%ind = spac%ihole(3)
        spac%cfg(1)%ihole(3)%spin = 'a'

        spac%coeff(2) = -0.5d0
        spac%cfg(2)%ipart(1)%ind = spac%ipart(1)
        spac%cfg(2)%ipart(1)%spin = 'b'
        spac%cfg(2)%ihole(1)%ind = spac%ihole(1)
        spac%cfg(2)%ihole(1)%spin = 'a'
        spac%cfg(2)%ihole(2)%ind = spac%ihole(2)
        spac%cfg(2)%ihole(2)%spin = 'b'
        spac%cfg(2)%ihole(3)%ind = spac%ihole(3)
        spac%cfg(2)%ihole(3)%spin = 'b'

        spac%coeff(3) = -0.5d0
        spac%cfg(3)%ipart(1)%ind = spac%ipart(1)
        spac%cfg(3)%ipart(1)%spin = 'a'
        spac%cfg(3)%ihole(1)%ind = spac%ihole(1)
        spac%cfg(3)%ihole(1)%spin = 'a'
        spac%cfg(3)%ihole(2)%ind = spac%ihole(2)
        spac%cfg(3)%ihole(2)%spin = 'a'
        spac%cfg(3)%ihole(3)%ind = spac%ihole(3)
        spac%cfg(3)%ihole(3)%spin = 'b'

        spac%coeff(4) = 0.5d0
        spac%cfg(4)%ipart(1)%ind = spac%ipart(1)
        spac%cfg(4)%ipart(1)%spin = 'b'
        spac%cfg(4)%ihole(1)%ind = spac%ihole(1)
        spac%cfg(4)%ihole(1)%spin = 'b'
        spac%cfg(4)%ihole(2)%ind = spac%ihole(2)
        spac%cfg(4)%ihole(2)%spin = 'b'
        spac%cfg(4)%ihole(3)%ind = spac%ihole(3)
        spac%cfg(4)%ihole(3)%spin = 'a'

      elseif(spac%typ==22) then
      endif

  endif

elseif(spac%spin==3) then

endif



end subroutine def_spac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module def_config

