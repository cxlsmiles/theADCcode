subroutine master_tda()

  use constants
  use parameters
  use sym_allowed
  use fspace
  use get_moment
  use misc
  use photoionisation
  
  implicit none

  integer, dimension(:,:), allocatable :: kpq
  integer :: i,ndim,nout 
  
  integer, dimension(norder) :: indx
  real(d), dimension(:), allocatable :: ener,mtm,tmvec,osc_str
  real(d), dimension(:,:), allocatable :: arr,domstates 

  allocate(kpq(7,0:nBas**2*4*nOcc**2))
  kpq(:,:)=-1

  call get_symallowed_tda(kpq(:,:))
  
  ndim=kpq(1,0)
  nout=ndim
!!$  if (ndim .gt. davstates) nout=davstates

  write(6,*) 'Direct diagonalisation of the TDA matrix'  
  
  allocate(arr(ndim,ndim),ener(ndim),mtm(ndim),tmvec(nout),osc_str(nout))
  call get_fspace_tda_direct(ndim,kpq(:,:),arr(:,:),ener(:))

  if (tranflag .eq. 'y') then
     write(6,*) 'Calculating transition moments in ',tranmom,' direction.'
     call get_modifiedtm_tda(ndim,kpq(:,:),mtm(:))
     do i= 1,nout
        tmvec(i)=tm(ndim,arr(:,i),mtm(:))
        osc_str(i)=2._d/3._d*ener(i)*tmvec(i)**2
        write(6,*) i,ener(i), os2cs*osc_str(i)
     end do
     
     call get_dominant(ndim,ener(:),osc_str(:),indx(:))
     allocate(domstates(ndim,norder))
     do i= 1,norder
        domstates(:,i)=arr(:,indx(i))
     end do
     call get_lanc_conf(ndim,domstates(:,:))
     deallocate(domstates)
     
!!$     call get_sums(ndim,ener(:),osc_str(:))
  end if

  if (tranflag .eq. 'y') then
     if (ndim .gt. davstates) nout=davstates
     call table2(ndim,nout,ener(1:nout),arr(:,1:nout),tmvec(:),osc_str(:))
  else
     if (ndim .gt. davstates) nout=davstates
     call table1(ndim,nout,ener(1:nout),arr(:,1:nout))
  end if

  call get_sigma(ndim,ener(:),os2cs*osc_str(:))  
  
!!  call get_symallowed_adc(kpq(:,:))

  deallocate(arr,ener,mtm,tmvec,osc_str,kpq)

end subroutine master_tda
     
  
  
