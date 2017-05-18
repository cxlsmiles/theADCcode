module fspace
  
  use constants
  use parameters
  use get_matrix
  use get_moment
  use select_fano
  use filetools
  use misc
  
  implicit none
  
contains

!!$----------------------------------------------------------

  subroutine get_fspace_adc2e_direct(ndim,kpq,arr,evector) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: ndim1,ndim2,nbuf,i,j,k
    
    real(d), dimension(:), allocatable :: ar_diag,temp
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    
    call get_offdiag_adc2ext_direct(ndim,kpq(:,:),ar_offdiag(:,:))
    call get_diag_adc2ext_direct(ndim1,ndim2,kpq(:,:),ar_diag(:))
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)

    call vdiagonalise(ndim,arr(:,:),evector(:))

  end subroutine get_fspace_adc2e_direct

!!$----------------------------------------------------------------------
!!$----------------------------------------------------------------------

  subroutine write_fspace_adc2e(ndim,kpq,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer*4 :: noffdel
    integer :: ndim1,ndim2,nbuf,i
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2ext_save(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2ext_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_adc2e

!!$-----------------------------------------------------------------------------------------
!!$----------------------------------------------------------------------

  subroutine write_fspace_adc2e_1(ndim,kpq,noffdel,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    integer*4, intent(out) :: noffdel

    integer :: ndim1,ndim2,nbuf,i
    logical :: prsd, prso
    
!!$    inquire(file='hmlt.dia'//chr,exist=prsd)
!!$    inquire(file='hmlt.off'//chr,exist=prso)
!!$    
!!$    if(prsd .and. prso) then
!!$       write(6,*) 'Engaging older hamiltonian files'
!!$       
!!$    else
    
       ndim1=kpq(1,0)
       ndim2=ndim-kpq(1,0)
       write(6,*) 'in write_space'
       call get_offdiag_adc2ext_save(ndim,kpq(:,:),nbuf,noffdel,chr)
       call get_diag_adc2ext_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
!!$    end if

  end subroutine write_fspace_adc2e_1
!!$-----------------------------------------------------------------------------------------

  subroutine get_fspace_adc2_direct(ndim,kpq,arr,evector) 
    
    integer :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: ndim1, ndim2, nbuf,i
    
    real(d), dimension(:), allocatable :: ar_diag
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    
    call get_offdiag_adc2_direct(ndim,kpq(:,:),ar_offdiag(:,:))
    call get_diag_adc2_direct(ndim1,ndim2,kpq(:,:),ar_diag(:))
    
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)
    
    call vdiagonalise(ndim,arr(:,:),evector(:))
    
    write(*,*) 'V Diag Done'

  end subroutine get_fspace_adc2_direct

!!$-----------------------------------------------------------------------------
!!$-----------------------------------------------------------------------------

  subroutine write_fspace_adc2(ndim,kpq,chr) 

    integer, intent(in) :: ndim 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer :: ndim1,ndim2,nbuf,i
    integer*4 :: noffdel
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2_save(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_adc2

!!$-----------------------------------------------------------------------------
!!$-----------------------------------------------------------------------------

  subroutine write_fspace_adc2_1(ndim,kpq,noffdel,chr) 

    integer, intent(in) :: ndim 
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr
    integer*4, intent(out) :: noffdel

    integer :: ndim1, ndim2, nbuf,i
    
    ndim1=kpq(1,0)
    ndim2=ndim-kpq(1,0)
    
    call get_offdiag_adc2_save(ndim,kpq(:,:),nbuf,noffdel,chr)
    call get_diag_adc2_save(ndim1,ndim2,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_adc2_1

!!$------------------------------------------------------------------------

!!$------------------------------------------------------------------------

  subroutine get_fspace_tda_direct(ndim,kpq,arr,evector) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq

    real(d), dimension(ndim), intent(out) :: evector
    real(d), dimension(ndim,ndim), intent(inout) :: arr

    integer :: nbuf,i
    
    real(d), dimension(:), allocatable :: ar_diag
    real(d), dimension(:,:), allocatable :: ar_offdiag
    
    allocate(ar_diag(ndim),ar_offdiag(ndim,ndim))
    call get_offdiag_tda_direct(ndim,kpq(:,:),ar_offdiag(:,:))
    call get_diag_tda_direct(ndim,kpq(:,:),ar_diag(:))
    
    
    arr(:,:)=ar_offdiag(:,:)
    
    do i=1,ndim
       arr(i,i)=ar_diag(i)
    end do
    
    deallocate(ar_diag,ar_offdiag)
    
    call vdiagonalise(ndim,arr(:,:),evector(:))
    
  end subroutine get_fspace_tda_direct

!!$-----------------------------------------------------------------------------
!!$-----------------------------------------------------------------------------

  subroutine write_fspace_tda(ndim,kpq,chr) 

    integer, intent(in) :: ndim
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(1), intent(in) :: chr

    integer :: nbuf,i
    
    call get_offdiag_tda_save(ndim,kpq(:,:),nbuf,chr)
    call get_diag_tda_save(ndim,kpq(:,:),nbuf,chr)
    
  end subroutine write_fspace_tda
  
!!$--------------------------------------------------------------------------
!!$--------------------------------------------------------------------------

  subroutine get_fstates_direct(ndim,kpq,fspace,fen,nstate,chflag1,chflag2)
    
    integer, intent(in) :: ndim
    integer, intent(out) :: nstate
    integer, dimension(7,0:nBas**2*nOcc**2), intent(in) :: kpq
    character(3), intent(in) :: chflag1
    character(1), intent(in) :: chflag2
    real(d), dimension(ndim,ndim), intent(out) :: fspace
    real(d), dimension(ndim), intent(out) :: fen
    
    integer :: nryd,ifail
    integer, dimension(nBas) :: ncnfi_ryd
    
    integer, dimension(:), allocatable :: nisri
    real(d), dimension(:,:), allocatable :: arr
    real(d), dimension(:), allocatable :: evec,temp
    
    integer :: i,lim
    integer*8 :: nstatel
    integer*8, dimension(:), allocatable :: nisril
    
    external M01CBF

    allocate(evec(ndim),arr(ndim,ndim),nisri(ndim),temp(ndim))
    if (chflag1 .eq. 'tda') then
       call get_fspace_tda_direct(ndim,kpq(:,:),arr(:,:),evec(:))
    elseif (chflag1 .eq. 'ad2') then
       call get_fspace_adc2_direct(ndim,kpq(:,:),arr(:,:),evec(:))
    elseif (chflag1 .eq. 'a2e') then
       call get_fspace_adc2e_direct(ndim,kpq(:,:),arr(:,:),evec(:))
    end if
    
!    if (chflag2 .eq. 'i') then
!       call get_ncnfi_ryd(kpq(:,:),ncnfi_ryd(:),nryd)
!       lim=nryd
!       temp(:)=0._d
!       do i=1,lim
!          temp(:)=temp(:)+arr(ncnfi_ryd(i),:)**2
!       end do
!    elseif (chflag2 .eq. 'f') then
!       lim=kpq(1,0)
!       temp(:)=0._d
!       do i=1,lim
!          temp(:)=temp(:)+arr(i,:)**2
!       end do
!    end if
    
!    if (chflag2 .eq. 'i') then
!       call select_fstate_ryd(ndim,mspacewi,temp(:),nstate,nisri(:))
!       write(6,*) nstate, "initial states has been selected"
!    elseif (chflag2 .eq. 'f') then
!       call select_fstate_ryd(ndim,mspacewf,temp(:),nstate,nisri(:))
!       write(6,*) nstate, "final states has been selected"
!    end if

     nstate=0

!    ifail=0
!    nstatel=int(nstate,lng)
!    nisril=int(nisri,lng)
!    call M01CBF(nisri(1:nstate),1,nstate,'A',ifail)

!    if ((chflag2 .eq. 'i') .and. (numinista .gt. nstate) ) then
!       write(6,*) "Requested number of the initial states is larger than the number of the available Rydberg states"
!       stop
!    end if

!    do i=1,nstate
!          write(6,*) nisri(i), evec(nisri(i))
!          fspace(:,i)=arr(:,nisri(i))
!          fen(i)=evec(nisri(i))
!       end do

    do i=1,ndim
          write(6,*) i, evec(i)
          fspace(:,i)=arr(:,i)
          fen(i)=evec(i)
       end do

    
    deallocate(evec,arr,nisri,temp)
    
  end subroutine get_fstates_direct
!!$--------------------------------------------------------------

  subroutine load_fstates_weight(ndim,negvc,kpq,fspace,fen,nstate,name,chflag)

    integer, intent(in) :: ndim,negvc
    integer, intent(out):: nstate
    integer, dimension(7,0:nBas**2*nOcc**2),intent(in) :: kpq 
    real(d), dimension(negvc), intent(out) :: fen
    real(d), dimension(ndim,negvc), intent(out) :: fspace
    character(36), intent(in) :: name
    character(1), intent(in) :: chflag

    integer, dimension(nBas) :: ncnfi_ryd
    integer, dimension(:), allocatable :: nisri,isv,indx
    real(d), dimension(:), allocatable :: ener,temp,ener_ryd
    real(d), dimension(:,:), allocatable :: arr

    integer :: i,vectype,vecdim,nvec,nryd,lim

    allocate(ener(negvc),arr(ndim,negvc),nisri(negvc),temp(negvc),isv(negvc))
    
    if (chflag .eq. 'i') then
       call readvct(ndim,1,negvc,ener(:),arr(:,:),nvec)
    elseif (chflag .eq. 'f') then
       call readvct(ndim,2,negvc,ener(:),arr(:,:),nvec)
    end if

    if(nvec .ne. negvc) then
       write(6,*) 'The number of read vectors',nvec,'differs from the number of requested states',negvc
       stop
    end if

    if (chflag .eq. 'i') then
       call get_ncnfi_ryd(kpq(:,:),ncnfi_ryd(:),nryd)
       lim=nryd
       temp(:)=0._d       
       do i=1,lim
          temp(:)=temp(:)+arr(ncnfi_ryd(i),:)**2
       end do
    elseif (chflag .eq. 'f') then
       lim=kpq(1,0)
       temp(:)=0._d
       do i=1,lim
          temp(:)=temp(:)+arr(i,:)**2
       end do
    end if

    if (chflag .eq. 'i') then
       call select_fstate_ryd(negvc,mspacewi,temp(:),nstate,nisri(:))
    elseif (chflag .eq. 'f') then
       call select_fstate_ryd(negvc,mspacewf,temp(:),nstate,nisri(:))
    end if
    write(6,*) nstate, "states has been selected"

    
    allocate(ener_ryd(nstate),indx(nstate))
    
    do i=1,nstate
       ener_ryd(i)=ener(nisri(i))
    end do
    
    call dsortindxa1('A',nstate,ener_ryd(:),indx(:))
    
    if (chflag .eq. 'i') then
       if (numinista .gt. nstate) then 
          write(6,*) "Number of the set initial states is larger than the number of the available Rydberg states"
          stop
       end if
    
       do i= 1,numinista
          fspace(:,i)=arr(:,nisri(indx(i)))
          fen(i)=ener(nisri(indx(i)))
       end do
       
    end if

    do i= 1,nstate
       fspace(:,i)=arr(:,nisri(indx(i)))
       fen(i)=ener(nisri(indx(i)))
    end do

    deallocate(ener,arr,nisri,temp,ener_ryd,indx)

  end subroutine load_fstates_weight
  
!!$--------------------------------------------------------------------------------
  
  subroutine get_bound(ndim,negvc,kpq,fspace,fen,nstate)
    
    integer, intent(in) :: ndim,negvc
    integer, intent(out):: nstate
    integer, dimension(7,0:nBas**2*nOcc**2),intent(in) :: kpq 
    real(d), dimension(negvc), intent(out) :: fen
    real(d), dimension(ndim,negvc), intent(out) :: fspace

    integer, dimension(nBas) :: ncnfi_ryd
    integer, dimension(:), allocatable :: nisri,isv,indx
    real(d), dimension(:), allocatable :: ener,temp,ener_ryd
    real(d), dimension(:,:), allocatable :: arr

    integer :: i,vectype,vecdim,nvec,nryd,lim

    allocate(ener(negvc),arr(ndim,negvc),nisri(negvc),temp(negvc),isv(negvc))
    
    !Reading ALL Davidson vectors 
    call readvct(ndim,1,negvc,ener(:),arr(:,:),nvec)

    call get_ncnfi_ryd(kpq(:,:),ncnfi_ryd(:),nryd)
    lim=nryd
    temp(:)=0._d       
    do i=1,lim
       temp(:)=temp(:)+arr(ncnfi_ryd(i),:)**2
    end do

    call select_fstate_ryd(negvc,mspacewi,temp(:),nstate,nisri(:))

    write(6,*) nstate, "states has been selected"

    
    allocate(ener_ryd(nstate),indx(nstate))
    
    do i=1,nstate
       ener_ryd(i)=ener(nisri(i))
    end do
    
    
    call dsortindxa1('A',nstate,ener_ryd(:),indx(:))
    
    if (numinista .gt. nstate) then 
       write(6,*) "Number of the set initial states is larger than the number of the available Rydberg states"
       stop
    end if
    
    do i= 1,numinista
       fspace(:,i)=arr(:,nisri(indx(i)))
       fen(i)=ener(nisri(indx(i)))
    end do
       
    do i= 1,nstate
       fspace(:,i)=arr(:,nisri(indx(i)))
       fen(i)=ener(nisri(indx(i)))
    end do

    deallocate(ener,arr,nisri,temp,ener_ryd,indx)

  end subroutine get_bound
    
!!$--------------------------------------------------------------

  subroutine get_tranmom_1(ndim,negvc,name,mtm,nstate,fen,tmvec)
    
    integer, intent(in) :: ndim,negvc
    integer, intent(out) :: nstate
    real(d), dimension(ndim), intent(in) :: mtm 
    real(d), dimension(negvc), intent(out) :: fen,tmvec
    
    character(36), intent(in) :: name

    integer :: i,num
    real(d) :: enr
    real(d), dimension(:), allocatable :: vec 
    logical :: log1

    INQUIRE(file=name,exist=log1)
    
    if(.not. log1) then
       write(6,*) 'The file ', name,' does not exist'
       stop
    end if
    
    allocate(vec(ndim))
    nstate=0
    
    OPEN(unit=77,file=name,status='OLD',access='SEQUENTIAL',form='UNFORMATTED')
    do i=1,negvc
       read(77,END=77) num,enr,vec(:)
       nstate=nstate+1
       write(6,*) 'Read vector',num,enr
       fen(nstate)=enr
       tmvec(nstate)=tm(ndim,vec(:),mtm(:))
    end do
77  CLOSE(77)
    write(6,*) 'There are ',nstate,' energies and tran. moments'
    
    
       
  end subroutine get_tranmom_1

!!$---------------------------------------------------------

!!$  subroutine get_gamma(ndim1,ndim2,fsten,gammavec)
!!$    
!!$    integer, intent(in) :: ndim1,ndim2
!!$    real(d), dimension(ndim1), intent(in) :: fsten
!!$    real(d), dimension(ndim1,ndim2), intent(in) :: gammavec
!!$
!!$    integer :: unitout,i,j
!!$    real(d) :: gamma0
!!$    character :: name
!!$    real(d), dimension(ndim1) :: dE
!!$
!!$    write(6,*) 'Running Stiltjes'
!!$
!!$    do i= 1,ndim2
!!$       unitout=10000+i
!!$           
!!$       do j= 1,ndim1
!!$          dE(j)=fsten(j)
!!$       end do
!!$           
!!$       call stieltjes (ndim1,dE(:),gammavec(:,i),gamma0,unitout,stiprilev)
!!$       name="state_"//achar(48+i)
!!$       open(10,file=name,access='sequential',status='unknown')
!!$       write(10,*) gamma0
!!$       close(10)
!!$           
!!$    end do
!!$
!!$  end subroutine get_gamma

!!$----------------------------------------------
  subroutine get_gamma_1(ndim,fsen,gammavec)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: fsen
    real(d), dimension(ndim), intent(in) :: gammavec

    integer :: i,ncount
    real(d) :: oslimit
    real(d), dimension(:), allocatable :: sgmvc_short,ener_short

    write(6,*) 'Printing Partial Gammas'
    allocate(sgmvc_short(ndim),ener_short(ndim))
    oslimit=1.e-16
    ncount=0
       
       do i= 1,ndim
          if (gammavec(i) .gt. oslimit) then
             ncount=ncount+1
             sgmvc_short(ncount)=gammavec(i)
             ener_short(ncount)=fsen(i)
          end if
       end do
       
       write(6,*) ncount,' states were estimated'
       do i=1,ncount
        write(6,*) ener_short(i),sgmvc_short(i)
       end do 
       
    deallocate(sgmvc_short,ener_short)
           
  end subroutine get_gamma_1
!!$---------------------------------------------------

  subroutine get_sigma(ndim,ener,sigmavec)

    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: ener
    real(d), dimension(ndim), intent(in) :: sigmavec

    integer :: i,nlimit,ncount
    real(d) :: gamma0,ehole,oslimit
    real(d), dimension(:), allocatable :: sgmvc_short,ener_short


    allocate(sgmvc_short(ndim),ener_short(ndim))
    write(6,*) 'Running Stiltjes'
    ehole=e(hinit)
    
    oslimit=1.e-7
    ncount=0
    do i= 1,ndim
       if (sigmavec(i) .gt. oslimit) then
          ncount=ncount+1
          sgmvc_short(ncount)=sigmavec(i)
          ener_short(ncount)=ener(i)
       end if
    end do
    call get_sums(ncount,ener_short(1:ncount),sgmvc_short(1:ncount)/os2cs)
  
       write(6,*) ncount,' states are to be sent to ST subroutine'
       do i=1,ncount
        write(6,*) ener_short(i),sgmvc_short(i)
       end do
!       write(6,*) 'Gamma for IS  Energy:',isen(ninista)
       nlimit=3000


    nlimit=3000

    if (ncount .ge. nlimit) then
!       call stieltjes_phi(ehole,omega,nlimit,ener_short(1:nlimit),sgmvc_short(1:nlimit),gamma0,7)
    else
!       call stieltjes_phi(ehole,omega,ncount,ener_short(1:ncount),sgmvc_short(1:ncount),gamma0,7)
    end if

    deallocate(sgmvc_short,ener_short)

  end subroutine get_sigma

!!$--------------------------------------------------------

  subroutine get_sums(ndim,ener,fosc)
    
    integer, intent(in) :: ndim
    real(d), dimension(ndim), intent(in) :: ener, fosc
    
    real(d), dimension(0:50):: sums
    integer :: i,j
    real(d) :: elev,flev,ratio

    sums(0:10)=0._d
    do i=1,ndim
       elev=ener(i)
       flev=fosc(i)
       ratio=flev
       do j=0,50
          sums(j)=sums(j)+ratio
          ratio=ratio/elev
       end do
    end do
    
    write(6,*) 'calculating the sums'
    write(6,*) 'S0 to S-50'
    do i= 0,50
       write(6,*) i,sums(i)
    end do
   
  end subroutine get_sums
!!$---------------------------------------------------------

  subroutine read_gammavec(ndimf,old_ls,nstatef,hpdimf,mgammavec,fsen,nosstates,nfsen,gammavec)
    
    integer, intent(in) :: ndimf,nstatef,old_ls
    integer, intent(out) :: nosstates
    real(d), dimension(old_ls), intent(in) :: fsen
    real(d), dimension(nstatef), intent(out) :: nfsen
    real(d), dimension(ndimf),intent(in) :: mgammavec
    real(d), dimension(nstatef), intent(out) :: gammavec
 
    integer :: mmr,nlim1,nlim2,fnm,lc1,lc2,nvecout,i,hpdimf
    integer :: lc1cnt
    real(d) :: pgamma,cntr,etest
    real(d), dimension(ndimf) :: tempvec1,tempvec2
    real(d), dimension(:,:), allocatable :: rvec 

    mmr=32768000!corr. to 250 Mb mem used in reading Lan. vecs from the saved file
    nlim1=mmr/ndimf
    nlim2=lancstates/nlim1+1
    fnm=2

    cntr=0.0
    lc1cnt=0
    nvecout=0
    allocate(rvec(ndimf,nlim1))
       tempvec1(:)=mgammavec(:)
       do lc2=1,nlim2 !counts the batches of read final states

          call readvct1(ndimf,fnm,(lc2-1)*nlim1+1,lc2*nlim1,rvec,nvecout)

          do lc1=1,nvecout !counts vectors inside the current batch

             do i=1,hpdimf
               cntr=cntr+rvec(i,lc1)**2
             end do
             etest=fsen((lc2-1)*nlim1+lc1)

!              if(etest.lt.0.5) then            ! 2h2p Variante : alle bis 0.5 mitnehmen (bis 0.5 wegen abs. anzahl)
!              write(6,*) 'etest',etest
!              lc1cnt=lc1cnt+1
!              tvecs(:,lc1cnt)=rvec(:,lc1)
!              tempvec2(:)=rvec(:,lc1)
!              pgamma=dsp(ndimf,tempvec1(:),tempvec2(:))
!              gammavec(lc1cnt,kc)=2._d*pi*pgamma**2        ! geht nur für lc1cnt << nlim1
!              nfsen(lc1cnt)=fsen((lc2-1)*nlim1+lc1)
!              else
!             write(6,*) 'etest',etest

             if(etest.lt.4) then
             if(cntr.gt.0.05) then
             lc1cnt=lc1cnt+1
             tempvec2(:)=rvec(:,lc1)
             pgamma=dsp(ndimf,tempvec1(:),tempvec2(:))
             gammavec(lc1cnt)=2._d*pi*pgamma**2
             nfsen(lc1cnt)=fsen((lc2-1)*nlim1+lc1)
             end if
             end if

             cntr=0.0
          end do
       end do

    nosstates=lc1cnt
    write(6,*) 'vectors read', (lc2-2)*nlim1+nvecout
    write(6,*) 'vectors taken', lc1cnt
    deallocate(rvec) 

  end subroutine read_gammavec
   
!!$---------------------------------------------------

  subroutine fill_stvc(ndim,vctr)
    
    integer, intent(in) :: ndim
    real(d),dimension(ndim), intent(in) :: vctr
    
    integer :: i,cnt
    integer, dimension(:), allocatable :: indarr
    
    allocate(indarr(ndim))
    call dsortindxa1('D',ndim,vctr(:)**2,indarr(:))

    stvc_lbl(1:lmain)=indarr(1:lmain)
    write(6,*)'Following configurations were chosen as the starting block of Lanzcos:',stvc_lbl(1:lmain)
    
    do i=1,lmain
       write(6,*) vctr(stvc_lbl(i))**2
    end do
    
    deallocate(indarr)
    
  end subroutine fill_stvc

!----------------------------------------------------------------------

 subroutine test_ortho(ndim,nstate,arr)

    integer, intent(in) :: ndim,nstate
    real(d), dimension(nstate,ndim), intent(in) :: arr
    real(d), dimension(:,:), allocatable :: mat
    real(d), dimension(:), allocatable :: tempvec1,tempvec2,evec
    integer :: i,j
    real(d) :: entri

    write(6,*) 'Starting Ortho Test'

    allocate(mat(ndim,ndim))
    allocate(tempvec1(nstate))
    allocate(tempvec2(nstate))

    do i=1,ndim
     do j=1,ndim
      tempvec1(:)=arr(:,i)
      tempvec2(:)=arr(:,j)
      entri=dsp(nstate,tempvec1(:),tempvec2(:))
      mat(i,j)=entri
     end do
    end do

    deallocate(tempvec1,tempvec2)

    do i=1,ndim
     write(*,'(99(f3.2,0X))') (mat(i,j),j=1,ndim)
    end do

    allocate(evec(ndim))

    call vdiagonalise(ndim,mat,evec)

    do i=1,ndim
     write(*,'(A,I3,2X,f4.3)') 'OM.EVal',i,evec(i)
    end do

    deallocate(mat,evec)

 end subroutine test_ortho

 subroutine show_vecs(ndim,nstate,arr)
! ndim=Anzahl Vecs, nstate=Länge Vecs, arr=Vecs
    integer, intent(in) :: ndim,nstate
    real(d), dimension(nstate,ndim), intent(in) :: arr
    integer :: i,j

    write(*,*) 'Starting Show_Vecs'

     do j=1,nstate
      write(*,'(I5,99(1X,f6.3))') j,(arr(j,i),i=1,15)
     end do

 end subroutine show_vecs

end module fspace
