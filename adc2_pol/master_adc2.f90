subroutine master_adc2()

  use constants
  use parameters
  use sym_allowed
  use davmod
  use fspace
  use get_moment
  use misc
  
  implicit none

  integer, dimension(:,:), allocatable :: kpq
  integer :: i,ndim,ndims,nout,nstates
  integer*4 :: noffd

  real(d) :: time
  real(d), dimension(:), allocatable :: ener,mtm,tmvec,osc_str
  real(d), dimension(:,:), allocatable :: arr 
  logical :: prs,prs1,prs2


  if (info .eq. 1) then
     allocate(kpq(7,0:nBas**2*4*nOcc**2))
     call get_symallowed_adc(kpq(:,:))
     deallocate(kpq)
     stop
  end if

  allocate(kpq(7,0:nBas**2*4*nOcc**2))
  kpq(:,:)=-1

  call get_symallowed_adc(kpq(:,:))

  ndims=kpq(1,0)
  ndim=kpq(1,0)+kpq(2,0)+kpq(3,0)+kpq(4,0)+2*kpq(5,0)
  nout=ndim
  
  if (chrun .eq. 'dire') then
     
     write(6,*) 'Direct diagonalisation of the ADC2 matrix'  
     

     allocate(arr(ndim,ndim),ener(ndim),mtm(ndim),tmvec(nout),osc_str(nout))
     call get_fspace_adc2_direct(ndim,kpq(:,:),arr(:,:),ener(:))

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
     
     deallocate(arr,ener)

  elseif (chrun .eq. 'save') then
     
     inquire(file='hmlt.diac',exist=prs)
     inquire(file='hmlt.offc',exist=prs1)
     inquire(file='fort.111',exist=prs2)
     
     if(prs .and. prs1 .and. prs2) then
        write(6,*) 'Older hmlt-file will be engaged'
        read(111,*) noffd
     else
        write(6,*) 'Saving complete ADC2 matrix in file'
        call  write_fspace_adc2_1(ndim,kpq(:,:),noffd,'c') 
        write(111,*) noffd
     end if

     call cpu_time(time)
     write(6,*) 'Time=',time," s"

     if (tranflag .eq. 'y') then
        allocate(mtm(ndim))
        inquire(file='mmnt',exist=prs)
        if (prs) then
           write(6,*) 'Older mmnt-file will be engaged'
           call read_vec(ndim,'mmnt',mtm(:))
        else
           write(6,*) 'Calculating transition moments in ',tranmom,' direction.'
           call get_modifiedtm_adc2(ndim,kpq(:,:),mtm(:))
           call write_vec(ndim,'mmnt',mtm(:))
        end if
     end if

     call cpu_time(time)
     write(6,*) 'Time=',time," s"
     write(6,*) 'check', mtm(1:ndims)
     
     deallocate(kpq)
        
     call fill_stvc(ndims,mtm(1:ndims))
     mtxidl='full'
     call master_lancdiag(ndim,noffd,'c')

     call cpu_time(time)
     write(6,*) 'Time=',time," s"
     
     if (tranflag .eq. 'y') then
        
        allocate(ener(lancstates),tmvec(lancstates))

!!$ ***lancstates is at most ncyclesXmain, usually we expect it to be lower than that***

        call get_tranmom_1(ndim,lancstates,lancname,mtm(:),nstates,ener(:),tmvec(:))
        
        allocate(osc_str(nstates))
        do i= 1,nstates
           osc_str(i)=2._d/3._d*ener(i)*tmvec(i)**2
        end do

        call cpu_time(time)
        write(6,*) 'Time=',time," s"

        call get_sigma(nstates,ener(1:nstates),os2cs*osc_str(:))
        deallocate(ener,tmvec,osc_str)
        
     else
        
        write(6,*) 'There is no alternative to the tranflag=y'

     end if
  end if

end subroutine master_adc2
