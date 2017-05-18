module davmod
  
  use constants
  use parameters
  
  implicit none 

  integer :: ndm,mem,main
  integer*4 :: noffdel
 
  real(d), dimension(:), allocatable :: offdiag,diag
  integer, dimension(:), allocatable :: indi, indip
  integer, dimension(:), allocatable :: stbl_index
  character(1) :: fl
  
contains
!!$---------------------------------------------------  
  subroutine master_dav(ndim,noffd,flag)
    
    integer, intent(in) :: ndim
    integer*4, intent(in) :: noffd
    character(1), intent(in) :: flag

    integer :: maxbl, nrec
    logical :: rstr,prs

    fl=flag
!!$    fl='i'
    noffdel=noffd
    ndm=ndim

    inquire(file='dav_vecs.'//mtxidd,exist=prs)
    if (prs) then
       write(6,*) 'Older file ','dav_vecs.'//mtxidd,' will be engaged'
    else
       allocate(offdiag(noffdel),indi(noffdel),indip(noffdel),diag(ndm))
       write(*,*) 'Alloc in Master Dav Successful'
       call read_matrix()
       write(*,*) 'Matrix read successfully'
!!$ The initial state correpsonds to the unit vec in the direction iv, i.e myb0=.false., main=1, i1,i2=1.
!!$ Parameter transp is always .true.

       call dnvini(rstr)
       if (rstr) then
          write(6,*) 'This is the Davidson restart run'
       else
          write(6,*) 'This is the Davidson startup run'
       endif
       
       mem=1000
       myb0=.true.
       transp=.true.
       call dinvop(ndm,dmain,mem,mtxidd,myb0,transp,rmtxhd,rmtxq1,rmtxhq1)
       
       deallocate(offdiag,indi,indip,diag)

    end if
  end subroutine master_dav

!!$-----------------------------------------------------------------

  subroutine master_lancdiag(ndim,noffd,flag)
    
    integer, intent(in) :: ndim
    integer*4, intent(in) :: noffd

    character(1), intent(in) :: flag
    integer :: i,j
    
    external blnczs

    fl=flag
    ndm=ndim
    main=lmain
    noffdel=noffd

    allocate(diag(ndm),offdiag(noffdel),indi(noffdel),indip(noffdel))
    call read_matrix()
 
    write(*,*) 'MMMMMMMMMMMMMM'
!    write(*,*) 'mtxidl',mtxidl
!    write(*,*) 'ndim',ndim
!    write(*,*) 'main',main
!    write(*,*) 'ncycles',ncycles
!    write(*,*) 'maxmem',maxmem
!    write(*,*) 'memx',memx
!    write(*,*) 'mode',mode
!    write(*,*) 'nprint',nprint
!    write(*,*) 'wthr',wthr
!    write(*,*) 'erange',erange(:)
!    write(*,*) 'unit',unit
!    write(*,*) 'iparm',iparm(:)
!    write(*,*) 'fparm',fparm(:)
!    write(*,*) '==============='     


    call blnczs(mtxidl,ndim,main,ncycles,maxmem,memx,mode,nprint,wthr,erange(:),unit,iparm(:),fparm(:),&
         mtxq1_l,rmtxhq1)
       
    deallocate(diag,offdiag,indi,indip)


  end subroutine master_lancdiag


!!$---------------------------------------------------

  subroutine  read_matrix()
    
    real(d) :: mtrl
    integer :: mark1,maxbl,nrec
    integer :: irec,nlim,i
    integer*4 :: count
    
    real(d), dimension(:), allocatable :: dgl,offdgl
    integer, dimension(:), allocatable :: oi,oj
    

!!$ Reading diagonal elements
    
    
    OPEN(UNIT=77,FILE='hmlt.dia'//fl,STATUS='OLD',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    
    read(77) maxbl,nrec
    read(77) diag(:)
    
    CLOSE(77)

    allocate(offdgl(maxbl),oi(maxbl),oj(maxbl))
    
    OPEN(UNIT=78,FILE='hmlt.off'//fl,STATUS='OLD',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')
    
    count=0
    do irec= 1,nrec
       read(78) offdgl(:),oi(:),oj(:),nlim
       offdiag(count+1:count+nlim)=offdgl(1:nlim)
       indi(count+1:count+nlim)=oi(1:nlim)
       indip(count+1:count+nlim)=oj(1:nlim)
       count=count+nlim
    end do
    close(78)
    deallocate(offdgl,oi,oj)
    
  end subroutine read_matrix


!!$-----------------------------------

   subroutine rmtxhd(dgl,work)

    real(d), dimension(ndm), intent(out) :: dgl
    real(d), dimension(mem), intent(inout) :: work

    OPEN(UNIT=12,FILE='hmlt.dia'//fl,STATUS='OLD',ACCESS='SEQUENTIAL',&
         FORM='UNFORMATTED')

    read(12) 
    read(12) dgl(:)

    CLOSE(12)
    
  end subroutine rmtxhd

!!$-----------------------------------------------------------------------

 subroutine rmtxq1(b,i1,i2,work)

   integer, intent(in) :: i1,i2
   real(d), dimension(ndm,i2-i1+1),intent(out) :: b
   real(d), dimension(mem),intent(inout) :: work

   integer :: i


   b(:,:)=rzero
   do i=i1,i2
      b(i,i-i1+1)=rone
   end do
   

  end subroutine rmtxq1
!$$---------------------------------------------------
  subroutine mtxq1_l(amx,bmx,i1,i2,work)

    implicit none

    real(d), dimension(ndm,i2-i1+1) :: bmx
    real(d) :: work,amx
    integer :: i1,i2,i

    bmx=0.0d0
   
    do i=i1,i2
       bmx(stvc_lbl(i),i-i1+1)=1.0d0
    end do

  end subroutine mtxq1_l
!!$------------------------------------------------

  subroutine  rmtxhq1(qv,zv,nq,work)

    implicit none

    integer, intent(in) :: nq
    real(d), dimension(nq,ndm) :: qv
    real(d), dimension(nq,ndm) :: zv
    real(d)  :: work

    real(d) :: mtrl
    integer :: mark1,ndim,maxbl,nrec,type
    integer*4 :: irec,nlim,i,j,irow,jcol

    zv=0.0d0

    do i= 1,ndm
       zv(:,i)=diag(i)*qv(:,i)
    end do
     
    do i= 1,noffdel
       irow=indi(i)
       jcol=indip(i)
       mtrl=offdiag(i)
       
       zv(:,irow)=zv(:,irow)+mtrl*qv(:,jcol)
       zv(:,jcol)=zv(:,jcol)+mtrl*qv(:,irow)
       
    end do
 
  end subroutine rmtxhq1

!!$-------------------------------------------

end module davmod
