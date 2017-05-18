subroutine master_adc2ext()

  use constants
  use parameters
  use sym_allowed
  use fspace
  use get_moment
  use misc
  
  implicit none

  integer, dimension(:,:), allocatable :: kpq
  integer :: i,ndim,nout,nstates,     j

  real(d), dimension(:), allocatable :: ener,mtm,tmvec,osc_str
  real(d), dimension(:,:), allocatable :: arr 


  real(d) :: summ,summ1,summ2

 if (info .eq. 1) then
     allocate(kpq(7,0:nBas**2*4*nOcc**2))
     call get_symallowed_adc(kpq(:,:))
     deallocate(kpq)
     stop
  end if
 
  allocate(kpq(7,0:nBas**2*4*nOcc**2))
  kpq(:,:)=-1

  call get_symallowed_adc(kpq(:,:))

  ndim=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)
  nout=ndim

  if (chrun .eq. 'dire') then
     
     write(6,*) 'Direct diagonalisation of the ADC2e matrix'  
     
     allocate(arr(ndim,ndim),ener(ndim),mtm(ndim),tmvec(nout),osc_str(nout))
     call get_fspace_adc2e_direct(ndim,kpq(:,:),arr(:,:),ener(:))

     if (tranflag .eq. 'y') then
        write(6,*) 'Calculating transition moments in ',tranmom,' direction.'
        call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:))

        do i= 1,nout
           tmvec(i)=tm(ndim,arr(:,i),mtm(:))
           osc_str(i)=2._d/3._d*ener(i)*tmvec(i)**2
        end do

     end if

     if (tranflag .eq. 'y') then
        if (ndim .gt. davstates) nout=davstates
        call table2(ndim,nout,ener(1:nout),arr(:,1:nout),tmvec(:),osc_str(:))
     else
        if (ndim .gt. davstates) nout=davstates
        call table1(ndim,nout,ener(1:nout),arr(:,1:nout))
     end if

     call get_sigma(ndim,ener(:),os2cs*osc_str(:))

     deallocate(arr,ener,mtm,tmvec,osc_str)
     
  elseif (chrun .eq. 'save') then 
     
     write(6,*) 'Saving complete ADC2e matrix in file'
     call  write_fspace_adc2e(ndim,kpq(:,:),'c')
     
     if (tranflag .eq. 'y') then
        allocate(mtm(ndim))
        write(6,*) 'Calculating transition moments in ',tranmom,' direction.'
        call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:))

        summ=0

        do j=1,ndim
           summ=summ+mtm(j)**2
        end do
        
        write(6,*) 'summ=',summ
        call write_vec(ndim,'mmnt',mtm(:))
        call write_vec(ndim,'stvc',mtm(:)/sqrt(summ))
     end if
     
  elseif (chrun .eq. 'read') then 
     
     allocate(arr(ndim,lancstates),ener(lancstates),mtm(ndim))
     
!!$     call load_fstates_band(ndim,lancstates,lancname,nstates,arr(:,:),ener(:))
     write(6,*) 'The following ',nstates,' final states were selected'
     
     allocate(tmvec(nstates),osc_str(nstates))

     if (tranflag .eq. 'y') then
        write(6,*) 'Reading transition moments in ',tranmom,' direction.'
        call read_vec(ndim,'mmnt',mtm(:))
        
!!$        write(6,*) 'Calculating transition moments in ',tranmom,' direction.'
!!$        call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:))
        
        do i= 1,nstates
           
           tmvec(i)=tm(ndim,arr(:,i),mtm(:))
           osc_str(i)=2._d/3._d*ener(i)*tmvec(i)**2
 
        end do

     end if

     nout=nstates
     if (nstates .gt. davstates) nout=davstates
     if (tranflag .eq. 'y') then
        call table2(ndim,nout,ener(1:nout),arr(:,1:nout),tmvec(:),osc_str(:))
     else
        call table1(ndim,nout,ener(1:nout),arr(:,1:nout))
     end if

     call get_sigma(nstates,ener(1:nstates),os2cs*osc_str(:))

     deallocate(arr,ener,mtm,tmvec,osc_str)
     
     
  end if  

end subroutine master_adc2ext
