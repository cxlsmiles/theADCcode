module read_param

  use constants
  use parameters
  use misc

  implicit none
  
  integer, parameter :: Bcknd=131072
  integer :: phis_init, ncheck
  
  external phis_init
  
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
subroutine read_molcas0()
  
  external phis_get_info
  
101 FORMAT(("|"),A25,3x,A25,3x,A25,("|"))
102 FORMAT(("|"),I25,3x,I25,3x,I25,("|"))
  
  ncheck=phis_init(Bcknd)

!!$ get nm. of irreps, basis functions, atomic centres
  
  call phis_get_info(nIrr,nBas,nCen)
  
  write(6,'(/,83("-"))')
  write(6,101) "Number of irreps","Number of active MO's", "Number of centres"
  write(6,102)  nIrr, nBas, nCen
  write(6,'(83("-")),/')
  

end subroutine read_molcas0
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_molcas1(info)  
  
  integer, intent(in) :: info
  integer :: j

  external phis_get_info
  external phis_get_epsi,phis_get_sym,phis_get_occ,phis_load_vpqrs

100 FORMAT(/,10x,A3,5x,A3,5x,A3,5x,A16)
101 FORMAT(/,10x,I3,5x,I3,5x,F3.1,5x,F16.10)
102 FORMAT(/,3("-"),A30,5x,F16.10,1x,A4)

  call phis_get_epsi(Ehf,e(:),nBas)
  call phis_get_sym(orbSym(:),nBas)
  call phis_get_occ(occNum(:),nBas)
  call phis_mc_dip(x_dipole(:,:),y_dipole(:,:),z_dipole(:,:),nBas)
  if(info.ne.1) then
  call phis_load_Vpqrs()
  end if
  write(6,102) "HF energy", Ehf, "a.u."
  write(6,100) 'nr.','sym','occ','orbital energy'
  write(6,'(10x,50("-"))')
  do j=1,nBas
     write(6,101) j,orbSym(j),occNum(j),e(j)
  end do
  

end subroutine read_molcas1

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine read_user()

  hcentre(:)=-1
  hneighb(:)=-1
  ninista=-1

  NAMELIST /USER/ nirrep,hcentre,hneighb,minc,stiprilev,method,idiag,fdiag,fmethod,&
       davname,lancname,mspacewi,mspacewf,davstates,numinista,chrun,&
       eupper,elower,readband,tranmom,norder,info,ninista
  NAMELIST /LNZLST/ ncycles,maxmem,memx,mode,nprint,maxiter,wthr,erange,unit,fparm,lmain

  READ(*,USER)
  READ(*,LNZLST)

  stvc_lbl(:)=-1

  iparm(1)=maxiter
  hinit=hcentre(1)
  lancstates=ncycles*lmain
  
  write(6,*) 'Orbitals on the central atom', hcentre(0),':',hcentre(1:hcentre(0))
  write(6,*) 'Orbitals on the neighb. atom', hneighb(0),':',hneighb(1:hneighb(0))
  
  write(6,*) "User data read"

end subroutine read_user


!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine check_user()
  

  if ((nirrep .lt. 1) .or. (nirrep .gt. 8)) then
     write(6,*) "IRREP IS NOT SET"
     stop
  elseif ((hinit .le. 0) .or. (hinit .gt. nBas)) then
     write(6,*) "INITIAL HOLE ORBITAL NUMBER IS NOT SET"
     stop
  elseif ((hinit .le. 0) .or. (hinit .gt. nBas)) then   
     write(6,*) "INITIAL PARTICLE ORBITAL NUMBER IS NOT SET"
     stop 
  elseif ((method .lt. 0) .or. (method .gt. 4)) then
     write(6,*) "COMPUTATION METHOD IS NOT CHOSEN"
     stop
  elseif ((numinista .lt. 0) .or. (numinista .gt. davstates)) then
     write(6,*) "THE NUMBER OF THE INITIAL STATES IS WRONG"
     stop   
  end if
     


end subroutine check_user

!------------------------------------------------------------------
!------------------------------------------------------------------

subroutine rearrange_occ()

!!$ The subroutine rearranges MO's such that occupied MO's precede unoccupied.
!!$ roccNum contains old indices, nOcc number of occup. orbs.
  
  integer :: i,j,k
  integer, dimension(nBas):: indx
  
100 FORMAT(/,3("-"),A50,3x,I4)
101 FORMAT(/,60("-"),/)

103 FORMAT(/,10x,A3,5x,A3,5x,A3,5x,A16)
104 FORMAT(/,10x,I3,5x,I3,5x,F3.1,5x,F16.10)
  
  j=0
  k=0
  do i=1,nBas
     if (nint(occNum(i)).eq.2) then
        j=j+1
!        roccNum(j)=i
     else
        k=k+1
!        roccNum(nAct-k+1)=i
     end if
  end do
  
  nOcc=j
  nVirt=k
  
  write(6,'(/,60("*"),/)')
  if ((k+j) .ne. nBas) then
     write(6,*) "(rearrange_occ):Orbital number is wrong. Bailing out ..."
     stop
  end if
  
  write(6,100) "hole space size",nOcc
  write(6,100) "particle space size",nVirt 
  write(6,101)
  
!!$  call dsortqx("A",nBas,e(:),1,indx(:))
  
  call dsortindxa1('A',nBas,e(:),indx(:))

  roccNum(:)=indx(:)
  
  write(6,103) 'nr.','sym','occ','orbital energy'
  write(6,'(10x,50("-"))')
  do i=1,nBas
     write(6,104) i,orbSym(roccNum(i)),occNum(roccNum(i)),e(roccNum(i))
  end do

end subroutine rearrange_occ
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module read_param
