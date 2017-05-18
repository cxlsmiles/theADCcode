module select_fano
  
  use constants
  use parameters
  use adc_ph
  use misc
  
  implicit none
  

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine select_singles(kpq,flag)
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in) :: flag
    
    integer :: i,ipart,ihole
    integer :: isym,a,j
    integer, dimension(7) :: col
    real(d) :: einit,ej
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))
    
    einit=abs(e(hinit))
    
    if (flag .eq. -1) then
       
       write(6,*) 'Selecting initial 1h1p subspace'
       
       kpq(1,0)=0
       
       do ihole=1,nOcc
          j=roccnum(ihole)
          ej=abs(e(j))
          do ipart=nOcc+1,nBas
             a=roccnum(ipart)
             isym=MT(orbSym(j),orbSym(a))
             if(isym .eq. nirrep) then
                if (einit .le. ej) then
                   kpq(1,0)=kpq(1,0)+1
                   call fill_indices(col(:),1,1,a,-1,j,-1,0) 
                   kpq(:,kpq(1,0))=col(:)
                end if
             end if
          end do
       end do
       
       write(6,100) "Number of 1h-1p configs in the IS ADC",kpq(1,0) 
       write(6,103)
       write(6,*) "Singly excited configurations allowed in the Init. St. Manif."
       write(6,101) "CNF","SPN","HL1","PT1"
       do i=1,kpq(1,0)
          write(6,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
       enddo
       
    elseif(flag .eq. 1) then
       
       write(6,*) 'Selecting final 1h1p subspace'
       
       kpq(1,0)=0
       
       do ihole=1,nOcc
          j=roccnum(ihole)
          ej=abs(e(j))
          do ipart=nOcc+1,nBas
             a=roccnum(ipart)
             isym=MT(orbSym(j),orbSym(a))
             if(isym .eq. nirrep) then
                if (einit .gt. ej) then
                   kpq(1,0)=kpq(1,0)+1
                   call fill_indices(col(:),1,1,a,-1,j,-1,0)
                   kpq(:,kpq(1,0))=col(:)
                end if
             end if
          end do
       end do
       
       write(6,100) "Number of 1h-1p configs in the FS ADC",kpq(1,0) 
       write(6,103)
       write(6,*) "Singly excited configurations allowed in the Fin. St. Manif."
       write(6,101) "CNF","SPN","HL1","PT1"
       do i=1,kpq(1,0)
          write(6,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
       enddo
       
    end if
    
  end subroutine select_singles
  
!!$-----------------------------------------------------------------------------------------------
  subroutine select_hetdia_inter_s(kpq)
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    
    integer :: i,ipart,ihole,n3smg
    integer :: isym,a
    integer, dimension(7) :: col
    
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))
    
    n3smg=hneighb(1)
    
    write(6,*) 'Selecting final 1h1p subspace for the intermolecular channel in hetdiatomic'
    
    kpq(1,0)=0
    
    do ipart=nOcc+1,nBas
       a=roccnum(ipart)
       isym=MT(orbSym(n3smg),orbSym(a))
       if(isym .eq. nirrep) then
          kpq(1,0)=kpq(1,0)+1
          call fill_indices(col(:),1,1,a,-1,n3smg,-1,0)
          kpq(:,kpq(1,0))=col(:)
       end if
    end do
       
    write(6,100) "Number of 1h-1p configs in the FS ADC",kpq(1,0) 
    write(6,103)
    write(6,*) "Singly excited configurations allowed in the Fin. St. Manif."
    write(6,101) "CNF","SPN","HL1","PT1"
    do i=1,kpq(1,0)
       write(6,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
    enddo
       
        
  end subroutine select_hetdia_inter_s
!!$-----------------------------------------------------------

  subroutine select_atom_is(kpq)
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    
    integer :: i,ap,n2pne
    integer :: isym,a,cah
    integer, dimension(7) :: col
    
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))
    
    write(6,*) 'Selecting initial 1h1p subspace'
    
    kpq(1,0)=0

    do ap=nOcc+1,nBas
       a=roccnum(ap)
       do i=1,hcentre(0)
           cah=hcentre(i)
          isym=MT(orbSym(cah),orbSym(a))
          if(isym .eq. nirrep) then
             kpq(1,0)=kpq(1,0)+1
             call fill_indices(col(:),1,1,a,-1,cah,-1,0)
             kpq(:,kpq(1,0))=col(:)
          end if
       end do
    end do
       
    write(6,100) "Number of 1h-1p configs in the IS ADC",kpq(1,0) 
    write(6,103)
    write(6,*) "Singly excited configurations allowed in the Fin. St. Manif."
    write(6,101) "CNF","SPN","HL1","PT1"
    do i=1,kpq(1,0)
       write(6,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
    enddo
       
        
  end subroutine select_atom_is

  subroutine select_atom_fs(kpq)
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    
    integer :: i,ap,n2pne
    integer :: isym,a,cah
    integer, dimension(7) :: col
    
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-"))
    
    write(6,*) 'Selecting final 1h1p subspace'
    
    kpq(1,0)=0

    do ap=nOcc+1,nBas
       a=roccnum(ap)
       do i=1,hneighb(0) 
           cah=hneighb(i)
          isym=MT(orbSym(cah),orbSym(a))
          if(isym .eq. nirrep) then
             kpq(1,0)=kpq(1,0)+1
             call fill_indices(col(:),1,1,a,-1,cah,-1,0)
             kpq(:,kpq(1,0))=col(:)
          end if
       end do
    end do
       
    write(6,100) "Number of 1h-1p configs in the IS ADC",kpq(1,0) 
    write(6,103)
    write(6,*) "Singly excited configurations allowed in the Fin. St. Manif."
    write(6,101) "CNF","SPN","HL1","PT1"
    do i=1,kpq(1,0)
       write(6,102) i,kpq(2,i),kpq(3,i),kpq(5,i)
    enddo
       
        
  end subroutine select_atom_fs

!!$-----------------------------------------------------------------------------------------------

    subroutine select_atom_d(kpq,flag)

!!$Includes 2p2pnxnx configurations in the both, initial and final states.
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in) :: flag

    integer :: i,j,a,b,k
    integer :: ih,jh,ap,bp
    integer :: cnti,cntf
    integer :: isym1, isym2
    integer, dimension(7) :: col
    real(d) :: einit,ei,ej

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 

    einit=abs(e(hinit))

! Filling the kpq arrays in the order: a=b,i=j;a|=b,i=j;a=b,i|=j;a|=b,i|=j(I,II).

    if (flag .eq. -1) then
       
       write(6,*) 'Selecting initial 2h2p subspace'
       
       kpq(2:5,0)=0
       cnti=kpq(1,0)
    
!!$a=b,i=j i=hcentre(ih)
       
       if(nirrep .eq. 1) then
          do ih=1,hcentre(0)
             i=hcentre(ih)
             ei=abs(e(i))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                if(einit .le. 2._d*ei) then
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                else
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                end if
             end do
          end do
       end if

!        if(nirrep .eq. 1) then
!           do ih=1,hneighb(0)
!              i=hneighb(ih)
!              ei=abs(e(i))
!              do ap=nOcc+1,nBas
!                 a=roccnum(ap)
!                 if(einit .le. 2._d*ei) then
!                    cnti=cnti+1
!                    kpq(2,0)=kpq(2,0)+1
!                    call fill_indices(col(:),2,1,a,a,i,i,1)
!                    kpq(:,cnti)=col(:)
!                 else
!                    cnti=cnti+1
!                    kpq(2,0)=kpq(2,0)+1
!                    call fill_indices(col(:),2,1,a,a,i,i,1)
!                    kpq(:,cnti)=col(:)
!                 end if
!              end do
!           end do
!        end if

!!$a|=b,i=j
  
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym1=MT(orbSym(a),orbSym(b))
             if(isym1 .eq. nirrep) then
                do ih=1,hcentre(0)
                   i=hcentre(ih)
                   ei=abs(e(i))
                   if(einit .le. 2._d*ei) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   else
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
          end do
       end do


!       do ap=nOcc+1,nBas
!          a=roccnum(ap)
!          do bp=ap+1,nBas
!             b=roccnum(bp)
!             isym1=MT(orbSym(a),orbSym(b))
!             if(isym1 .eq. nirrep) then
!                do ih=1,hneighb(0)
!                   i=hneighb(ih)
!                   ei=abs(e(i))
!                   if(einit .le. 2._d*ei) then
!                      cnti=cnti+1
!                      kpq(3,0)=kpq(3,0)+1
!                      call fill_indices(col(:),2,1,a,b,i,i,2)
!                      kpq(:,cnti)=col(:)
!                   else
!                      cnti=cnti+1
!                      kpq(3,0)=kpq(3,0)+1
!                      call fill_indices(col(:),2,1,a,b,i,i,2)
!                      kpq(:,cnti)=col(:)
!                   end if
!                end do
!             end if
!          end do
!       end do
!a=b,i|=j
  
       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             if (isym1 .eq. nirrep) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                end do
             end if
          end do
       end do

!       do ih=1,hneighb(0)
!          i=hneighb(ih)
!          ei=abs(e(i))
!          do jh=1,hcentre(0)
!             j=hcentre(jh)
!             ej=abs(e(j))
!             isym1=MT(orbSym(i),orbSym(j))
!             if (isym1 .eq. nirrep) then
!                do ap=nOcc+1,nBas
!                   a=roccnum(ap)
!                      cnti=cnti+1
!                      kpq(4,0)=kpq(4,0)+1
!                      call fill_indices(col(:),2,1,a,a,i,j,3) 
!                      kpq(:,cnti)=col(:)
!                end do
!             end if
!          end do
!       end do

!       do ih=1,hneighb(0)
!          i=hneighb(ih)
!          ei=abs(e(i))
!          do jh=ih+1,hneighb(0)
!             j=hneighb(jh)
!             ej=abs(e(j))
!             isym1=MT(orbSym(i),orbSym(j))
!             if (isym1 .eq. nirrep) then
!                do ap=nOcc+1,nBas
!                   a=roccnum(ap)
!                      cnti=cnti+1
!                      kpq(4,0)=kpq(4,0)+1
!                      call fill_indices(col(:),2,1,a,a,i,j,3) 
!                      kpq(:,cnti)=col(:)
!                end do
!             end if
!          end do
!       end do
    
!a|=b,i|=j spin I

       do ih=1,hcentre(0)
          i=hcentre(ih)
          ei=abs(e(i))
          do jh=ih+1,hcentre(0)
             j=hcentre(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                do bp=ap+1,nBas
                   b=roccnum(bp)
                   isym2=MT(orbSym(a),orbSym(b))
                   if(MT(isym1,isym2) .eq. nirrep) then 
                      if(einit .le. (ei+ej)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      else
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      end if
                   end if
                end do
             end do
          end do
       end do

!       do ih=1,hneighb(0)
!          i=hneighb(ih)
!          ei=abs(e(i))
!          do jh=1,hcentre(0)
!             j=hcentre(jh)
!             ej=abs(e(j))
!             isym1=MT(orbSym(i),orbSym(j))
!             do ap=nOcc+1,nBas
!                a=roccnum(ap)
!                do bp=ap+1,nBas
!                   b=roccnum(bp)
!                   isym2=MT(orbSym(a),orbSym(b))
!                   if(MT(isym1,isym2) .eq. nirrep) then 
!                      if(einit .le. (ei+ej)) then
!                         cnti=cnti+1
!                         kpq(5,0)=kpq(5,0)+1
!                         call fill_indices(col(:),2,11,a,b,i,j,4)
!                         kpq(:,cnti)=col(:)
!                      else
!                         cnti=cnti+1
!                         kpq(5,0)=kpq(5,0)+1
!                         call fill_indices(col(:),2,11,a,b,i,j,4)
!                         kpq(:,cnti)=col(:)
!                      end if
!                   end if
!                end do
!             end do
!          end do
!       end do

! aukommentier -> Ne sonst HF

!       do ih=1,hneighb(0)
!          i=hneighb(ih)
!          ei=abs(e(i))
!          do jh=ih+1,hneighb(0)
!             j=hneighb(jh)
!             ej=abs(e(j))
!             isym1=MT(orbSym(i),orbSym(j))
!             do ap=nOcc+1,nBas
!                a=roccnum(ap)
!                do bp=ap+1,nBas
!                   b=roccnum(bp)
!                   isym2=MT(orbSym(a),orbSym(b))
!                   if(MT(isym1,isym2) .eq. nirrep) then 
!                      if(einit .le. (ei+ej)) then
!                         cnti=cnti+1
!                         kpq(5,0)=kpq(5,0)+1
!                         call fill_indices(col(:),2,11,a,b,i,j,4)
!                         kpq(:,cnti)=col(:)
!                      else
!                         cnti=cnti+1
!                         kpq(5,0)=kpq(5,0)+1
!                         call fill_indices(col(:),2,11,a,b,i,j,4)
!                         kpq(:,cnti)=col(:)
!                      end if
!                   end if
!                end do
!             end do
!          end do
!       end do
    
       kpq(:,cnti+1:cnti+kpq(5,0))=kpq(:,cnti+1-kpq(5,0):cnti)
       kpq(2,cnti+1:cnti+kpq(5,0))=12
       kpq(7,cnti+1:cnti+kpq(5,0))=5

       
       write(6,100) "Number of 2h-2p |abij> configs in the IS ADC", cnti+kpq(5,0)-kpq(1,0)
       write(6,103)
       write(6,*) " Doubly excited |abij> configs allowed in the Init. St. Manif."
       write(6,101) "CNF","SPN","HL1","HL2","PT1","PT2"
    
       do k=kpq(1,0)+1,cnti+kpq(5,0) 
          write(6,102) k,kpq(2,k),kpq(3,k),&
               kpq(4,k),kpq(5,k),kpq(6,k)
       end do

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
    elseif(flag .eq. 1) then
       
       write(6,*) 'Selecting final 2h2p subspace'
       kpq(2:5,0)=0
       cntf=kpq(1,0)

!!$a=b,i=j
    
       if(nirrep .eq. 1) then
          do ih=1,hneighb(0)
             i=hneighb(ih)
             ei=abs(e(i))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
!                if(einit .gt. 2._d*ei) then
                   cntf=cntf+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cntf)=col(:)
!                end if
             end do
          end do
       end if

!       if(nirrep .eq. 1) then
!          do ih=3,4 !hcentre(0)
!             i=hcentre(ih)
!             ei=abs(e(i))
!             do ap=nOcc+1,nBas
!                a=roccnum(ap)
!                if(einit .gt. 2._d*ei) then
!                   cntf=cntf+1
!                   kpq(2,0)=kpq(2,0)+1
!                   call fill_indices(col(:),2,1,a,a,i,i,1)
!                   kpq(:,cntf)=col(:)
!                end if
!             end do
!          end do
!       end if

!a|=b,i=j
  

       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym1=MT(orbSym(a),orbSym(b))
             if(isym1 .eq. nirrep) then
                do ih=1,hneighb(0)
                   i=hneighb(ih)
                   ei=abs(e(i))
!                   if(einit .gt. 2._d*ei) then
                      cntf=cntf+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2) 
                      kpq(:,cntf)=col(:)
!                   end if
                end do
             end if
          end do
       end do

!       do ap=nOcc+1,nBas
!          a=roccnum(ap)
!          do bp=ap+1,nBas
!             b=roccnum(bp)
!             isym1=MT(orbSym(a),orbSym(b))
!             if(isym1 .eq. nirrep) then
!               do ih=1,hcentre(0)
!                   i=hcentre(ih)
!                   ei=abs(e(i))
!                   if(einit .gt. 2._d*ei) then
!                      cntf=cntf+1
!                      kpq(3,0)=kpq(3,0)+1
!                      call fill_indices(col(:),2,1,a,b,i,i,2) 
!                      kpq(:,cntf)=col(:)
!                   end if
!                end do
!             end if
!          end do
!       end do
       
!a=b,i|=j
  
       do ih=1,hneighb(0)
          i=hneighb(ih)
          ei=abs(e(i))
          do jh=ih+1,hneighb(0)
             j=hneighb(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             if (isym1 .eq. nirrep) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
!                   if(einit .gt. (ei+ej)) then
                      cntf=cntf+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cntf)=col(:)
!                   end if
                end do
             end if
          end do
       end do

!       do ih=1,hcentre(0)
!          i=hcentre(ih)
!          ei=abs(e(i))
!          do jh=1,hneighb(0)
!             j=hneighb(jh)
!             ej=abs(e(j))
!             isym1=MT(orbSym(i),orbSym(j))
!             if (isym1 .eq. nirrep) then
!                do ap=nOcc+1,nBas
!                   a=roccnum(ap)
!                   if(einit .gt. (ei+ej)) then
!                      cntf=cntf+1
!                      kpq(4,0)=kpq(4,0)+1
!                      call fill_indices(col(:),2,1,a,a,i,j,3) 
!                      kpq(:,cntf)=col(:)
!                   end if
!                end do
!             end if
!          end do
!       end do

!       do ih=1,hcentre(0)
!          i=hcentre(ih)
!          ei=abs(e(i))
!          do jh=ih+1,hcentre(0)
!             j=hcentre(jh)
!             ej=abs(e(j))
!             isym1=MT(orbSym(i),orbSym(j))
!             if (isym1 .eq. nirrep) then
!                do ap=nOcc+1,nBas
!                   a=roccnum(ap)
!                   if(einit .gt. (ei+ej)) then
!                      cntf=cntf+1
!                      kpq(4,0)=kpq(4,0)+1
!                      call fill_indices(col(:),2,1,a,a,i,j,3) 
!                      kpq(:,cntf)=col(:)
!                   end if
!                end do
!             end if
!          end do
!       end do
    
!a|=b,i|=j spin I

       do ih=1,hneighb(0)
          i=hneighb(ih)
          ei=abs(e(i))
          do jh=ih+1,hneighb(0)
             j=hneighb(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                do bp=ap+1,nBas
                   b=roccnum(bp)
                   isym2=MT(orbSym(a),orbSym(b))
                   if(MT(isym1,isym2) .eq. nirrep) then 
!                      if(einit .gt. (ei+ej)) then
                         cntf=cntf+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cntf)=col(:)
!                      end if
                   end if
                end do
             end do
          end do
       end do

!       do ih=1,hcentre(0)
!          i=hcentre(ih)
!          ei=abs(e(i))
!          do jh=1,hneighb(0)
!             j=hneighb(jh)
!             ej=abs(e(j))
!             isym1=MT(orbSym(i),orbSym(j))
!             do ap=nOcc+1,nBas
!                a=roccnum(ap)
!                do bp=ap+1,nBas
!                   b=roccnum(bp)
!                   isym2=MT(orbSym(a),orbSym(b))
!                   if(MT(isym1,isym2) .eq. nirrep) then 
!                      if(einit .gt. (ei+ej)) then
!                         cntf=cntf+1
!                         kpq(5,0)=kpq(5,0)+1
!                         call fill_indices(col(:),2,11,a,b,i,j,4)
!                         kpq(:,cntf)=col(:)
!                      end if
!                   end if
!                end do
!             end do
!          end do
!       end do

!       do ih=1,hcentre(0)
!          i=hcentre(ih)
!          ei=abs(e(i))
!         do jh=ih+1,hcentre(0)
!             j=hcentre(jh)
!             ej=abs(e(j))
!             isym1=MT(orbSym(i),orbSym(j))
!             do ap=nOcc+1,nBas
!                a=roccnum(ap)
!                do bp=ap+1,nBas
!                   b=roccnum(bp)
!                  isym2=MT(orbSym(a),orbSym(b))
!                   if(MT(isym1,isym2) .eq. nirrep) then 
!                      if(einit .gt. (ei+ej)) then
!                         cntf=cntf+1
!                         kpq(5,0)=kpq(5,0)+1
!                         call fill_indices(col(:),2,11,a,b,i,j,4)
!                         kpq(:,cntf)=col(:)
!                      end if
!                   end if
!                end do
!             end do
!          end do
!       end do
    
       kpq(:,cntf+1:cntf+kpq(5,0))=kpq(:,cntf+1-kpq(5,0):cntf)
       kpq(2,cntf+1:cntf+kpq(5,0))=12 
       kpq(7,cntf+1:cntf+kpq(5,0))=5
    
       write(6,100) "Number of 2h-2p |abij> configs in the FS ADC", cntf+kpq(5,0)-kpq(1,0)
       write(6,103)
       write(6,*) " Doubly excited |abij> configs allowed in the Fin. St. Manif."
       write(6,101) "CNF","SPN","HL1","HL2","PT1","PT2"
       
       do k=kpq(1,0)+1, cntf+kpq(5,0)
          write(6,102) k,kpq(2,k),kpq(3,k),&
               kpq(4,k),kpq(5,k),kpq(6,k)
       end do
       
    end if
    
  end subroutine select_atom_d
!!$--------------------------------------------


  subroutine select_doubles_total(kpq,flag)

!!$Includes 2p2pnxnx configurations in the both, initial and final states. Used for heteroatomic diatomics.
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    integer, intent(in) :: flag

    integer :: i,j,a,b,k
    integer :: ih,jh,ap,bp
    integer :: cnti,cntf,n3smg
    integer :: isym1, isym2
    integer, dimension(7) :: col
    real(d) :: einit,ei,ej,e2pne

100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 

    einit=abs(e(hinit))
    e2pne=abs(e(hcentre(2)))
    n3smg=hneighb(1)


! Filling the kpq arrays in the order: a=b,i=j;a|=b,i=j;a=b,i|=j;a|=b,i|=j(I,II).

    if (flag .eq. -1) then
       
       write(6,*) 'Selecting initial 1h1p subspace'
       
       kpq(2:5,0)=0
       cnti=kpq(1,0)
    
!!$a=b,i=j
       
       if(nirrep .eq. 1) then
          do ih=1,nOcc
             i=roccnum(ih)
             ei=abs(e(i))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                if(einit .le. 2._d*ei) then
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                elseif(i .ne. n3smg) then
                   cnti=cnti+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cnti)=col(:)
                end if
             end do
          end do
       end if

!!$a|=b,i=j
  
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym1=MT(orbSym(a),orbSym(b))
             if(isym1 .eq. nirrep) then
                do ih=1,nOcc
                   i=roccnum(ih)
                   ei=abs(e(i))
                   if(einit .le. 2._d*ei) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   elseif(i .ne. n3smg) then
                      cnti=cnti+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2)
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
          end do
       end do
    
!a=b,i|=j
  
       do ih=1,nOcc
          i=roccnum(ih)
          ei=abs(e(i))
          do jh=ih+1,nOcc
             j=roccnum(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             if (isym1 .eq. nirrep) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                   if((einit .le. (ei+ej)) .and. (i .ne. n3smg) .and. (j .ne. n3smg)) then
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                   elseif((i .ne. n3smg) .and. (j .ne. n3smg)) then
                      cnti=cnti+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cnti)=col(:)
                   end if
                end do
             end if
          end do
       end do
    
!a|=b,i|=j spin I

       do ih=1,nOcc
          i=roccnum(ih)
          ei=abs(e(i))
          do jh=ih+1,nOcc
             j=roccnum(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                do bp=ap+1,nBas
                   b=roccnum(bp)
                   isym2=MT(orbSym(a),orbSym(b))
                   if(MT(isym1,isym2) .eq. nirrep) then 
                      if((einit .le. (ei+ej)) .and. (i .ne. n3smg) .and. (j .ne. n3smg)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      elseif((i .ne. n3smg) .and. (j .ne. n3smg)) then
                         cnti=cnti+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cnti)=col(:)
                      end if
                   end if
                end do
             end do
          end do
       end do
    
       kpq(:,cnti+1:cnti+kpq(5,0))=kpq(:,cnti+1-kpq(5,0):cnti)
       kpq(2,cnti+1:cnti+kpq(5,0))=12
       kpq(7,cnti+1:cnti+kpq(5,0))=5

       
       write(6,100) "Number of 2h-2p |abij> configs in the IS ADC", cnti+kpq(5,0)-kpq(1,0)
       write(6,103)
       write(6,*) " Doubly excited |abij> configs allowed in the Init. St. Manif."
       write(6,101) "CNF","SPN","HL1","HL2","PT1","PT2"
    
       do k=kpq(1,0)+1,cnti+kpq(5,0) 
          write(6,102) k,kpq(2,k),kpq(3,k),&
               kpq(4,k),kpq(5,k),kpq(6,k)
       end do

    elseif(flag .eq. 1) then
       
       write(6,*) 'Selecting final 1h1p subspace'
       kpq(2:5,0)=0
       cntf=kpq(1,0)

!!$a=b,i=j
    
       if(nirrep .eq. 1) then
          do ih=1,nOcc
             i=roccnum(ih)
             ei=abs(e(i))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                if(einit .gt. 2._d*ei) then
                   cntf=cntf+1
                   kpq(2,0)=kpq(2,0)+1
                   call fill_indices(col(:),2,1,a,a,i,i,1)
                   kpq(:,cntf)=col(:)
                end if
             end do
          end do
       end if

!a|=b,i=j
  

       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym1=MT(orbSym(a),orbSym(b))
             if(isym1 .eq. nirrep) then
                do ih=1,nOcc
                   i=roccnum(ih)
                   ei=abs(e(i))
                   if(einit .gt. 2._d*ei) then
                      cntf=cntf+1
                      kpq(3,0)=kpq(3,0)+1
                      call fill_indices(col(:),2,1,a,b,i,i,2) 
                      kpq(:,cntf)=col(:)
                   end if
                end do
             end if
          end do
       end do
       
!a=b,i|=j
  
       do ih=1,nOcc
          i=roccnum(ih)
          ei=abs(e(i))
          do jh=ih+1,nOcc
             j=roccnum(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             if (isym1 .eq. nirrep) then
                do ap=nOcc+1,nBas
                   a=roccnum(ap)
                   if(einit .gt. (ei+ej)) then
                      cntf=cntf+1
                      kpq(4,0)=kpq(4,0)+1
                      call fill_indices(col(:),2,1,a,a,i,j,3) 
                      kpq(:,cntf)=col(:)
                   end if
                end do
             end if
          end do
       end do
    
!a|=b,i|=j spin I

       do ih=1,nOcc
          i=roccnum(ih)
          ei=abs(e(i))
          do jh=ih+1,nOcc
             j=roccnum(jh)
             ej=abs(e(j))
             isym1=MT(orbSym(i),orbSym(j))
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                do bp=ap+1,nBas
                   b=roccnum(bp)
                   isym2=MT(orbSym(a),orbSym(b))
                   if(MT(isym1,isym2) .eq. nirrep) then 
                      if(einit .gt. (ei+ej)) then
                         cntf=cntf+1
                         kpq(5,0)=kpq(5,0)+1
                         call fill_indices(col(:),2,11,a,b,i,j,4)
                         kpq(:,cntf)=col(:)
                      end if
                   end if
                end do
             end do
          end do
       end do
    
       kpq(:,cntf+1:cntf+kpq(5,0))=kpq(:,cntf+1-kpq(5,0):cntf)
       kpq(2,cntf+1:cntf+kpq(5,0))=12 
       kpq(7,cntf+1:cnti+kpq(5,0))=5
    
       write(6,100) "Number of 2h-2p |abij> configs in the FS ADC", cntf+kpq(5,0)-kpq(1,0)
       write(6,103)
       write(6,*) " Doubly excited |abij> configs allowed in the Fin. St. Manif."
       write(6,101) "CNF","SPN","HL1","HL2","PT1","PT2"
       
       do k=kpq(1,0)+1, cntf+kpq(5,0)
          write(6,102) k,kpq(2,k),kpq(3,k),&
               kpq(4,k),kpq(5,k),kpq(6,k)
       end do
       
    end if
    
  end subroutine select_doubles_total

!!$------------------------------------------------------

  subroutine select_hetdia_inter_d(kpq)
    
!!$Includes 2p2pnxnx configurations in the both, initial and final states. Used for heteroatomic diatomics.
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    
    integer :: i,j,a,b,k
    integer :: ih,jh,ap,bp
    integer :: cnti,cntf,n3smg,ncentre
    integer :: isym1, isym2
    integer, dimension(7) :: col
    
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 
    
    ncentre=hcentre(0)
    n3smg=hneighb(1)

!!$as the confs for the interatomic channel we choose: NeMg+ (in s-routine), Ne+Mg+, and NeMg++.
!!$Ne++Mg, which might be useful to describe CI in Ne*Mg+ are not taken into account.


! Filling the kpq arrays in the order: a=b,i=j;a|=b,i=j;a=b,i|=j;a|=b,i|=j(I,II).

       
    write(6,*) 'Selecting final 2h2p subspace for the intermolecular channel in hetdiatomic '
    kpq(2:5,0)=0
    cntf=kpq(1,0)
    
!!$a=b,i=j
    
    if(nirrep .eq. 1) then !!NeMg++,i=n3smg
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          cntf=cntf+1
          kpq(2,0)=kpq(2,0)+1
          call fill_indices(col(:),2,1,a,a,n3smg,n3smg,1)
          kpq(:,cntf)=col(:)
       end do
    end if
    
!a|=b,i=j  ; NeMg++,i=n3smg
    
    do ap=nOcc+1,nBas
       a=roccnum(ap)
       do bp=ap+1,nBas
          b=roccnum(bp)
          isym1=MT(orbSym(a),orbSym(b))
          if(isym1 .eq. nirrep) then
             cntf=cntf+1
             kpq(3,0)=kpq(3,0)+1
             call fill_indices(col(:),2,1,a,b,n3smg,n3smg,2) 
             kpq(:,cntf)=col(:)
          end if
       end do
    end do
       
!a=b,i|=j Ne+Mg+  i on Ne, j=n3sMg
  
    do ih=2,ncentre
       i=hcentre(ih)
       isym1=MT(orbSym(i),orbSym(n3smg))
       if (isym1 .eq. nirrep) then
          do ap=nOcc+1,nBas
             a=roccnum(ap)
             cntf=cntf+1
             kpq(4,0)=kpq(4,0)+1
             call fill_indices(col(:),2,1,a,a,i,n3smg,3) 
             kpq(:,cntf)=col(:)
          end do
       end if
    end do
    
!a|=b,i|=j spin I Ne+Mg+ i on Ne, j=n3sMg

    do ih=2,ncentre
       i=hcentre(ih)
       isym1=MT(orbSym(i),orbSym(n3smg))
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do bp=ap+1,nBas
             b=roccnum(bp)
             isym2=MT(orbSym(a),orbSym(b))
             if(MT(isym1,isym2) .eq. nirrep) then 
                cntf=cntf+1
                kpq(5,0)=kpq(5,0)+1
                call fill_indices(col(:),2,11,a,b,i,n3smg,4)
                kpq(:,cntf)=col(:)
             end if
          end do
       end do
    end do
    
    kpq(:,cntf+1:cntf+kpq(5,0))=kpq(:,cntf+1-kpq(5,0):cntf)
    kpq(2,cntf+1:cntf+kpq(5,0))=12 
    kpq(7,cntf+1:cnti+kpq(5,0))=5
    
    write(6,100) "Number of 2h-2p |abij> configs in the FS ADC", cntf+kpq(5,0)-kpq(1,0)
    write(6,103)
    write(6,*) " Doubly excited |abij> configs allowed in the Fin. St. Manif."
    write(6,101) "CNF","SPN","HL1","HL2","PT1","PT2"
    
    do k=kpq(1,0)+1, cntf+kpq(5,0)
       write(6,102) k,kpq(2,k),kpq(3,k),&
            kpq(4,k),kpq(5,k),kpq(6,k)
    end do
    
    
  end subroutine select_hetdia_inter_d

!!$------------------------------------------------------

  subroutine select_hetdia_au_d(kpq)
    
!!$Includes 2p2pnxnx configurations in the both, initial and final states. Used for heteroatomic diatomics.
    
    integer, dimension(7,0:nBas**2*nOcc**2), intent(inout) :: kpq
    
    integer :: i,j,a,b,k
    integer :: ih,jh,ap,bp
    integer :: cnti,cntf,n3smg,ncentre
    integer :: isym1, isym2
    integer, dimension(7) :: col
    
    
100 FORMAT(/,3("*"),A50,3x,I4)
101 FORMAT(/,("*"),3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,3x,A3,/)
102 FORMAT(("*"),3x,I5,3x,I3,3x,I3,3x,I3,3x,I3,3x,I3)
103 FORMAT(/,60("-")) 
    
    ncentre=hcentre(0)

!!$as the confs for the interatomic channel we choose: NeMg+ (in s-routine), Ne+Mg+, and NeMg++.
!!$Ne++Mg, which might be useful to describe CI in Ne*Mg+ are not taken into account.


! Filling the kpq arrays in the order: a=b,i=j;a|=b,i=j;a=b,i|=j;a|=b,i|=j(I,II).

       
    write(6,*) 'Selecting final 2h2p subspace for the intermolecular channel in hetdiatomic '
    kpq(2:5,0)=0
    cntf=kpq(1,0)
    
!!$a=b,i=j
    
    if(nirrep .eq. 1) then !!Ne++Mg
       do ap=nOcc+1,nBas
          a=roccnum(ap)
          do ih=2,ncentre
             i=hcentre(ih)
             cntf=cntf+1
             kpq(2,0)=kpq(2,0)+1
             call fill_indices(col(:),2,1,a,a,i,i,1)
             kpq(:,cntf)=col(:)
          end do
       end do
    end if
    
!a|=b,i=j 
    
    do ap=nOcc+1,nBas
       a=roccnum(ap)
       do bp=ap+1,nBas
          b=roccnum(bp)
          isym1=MT(orbSym(a),orbSym(b))
          if(isym1 .eq. nirrep) then
             do ih=2,ncentre
                i=hcentre(ih)
                cntf=cntf+1
                kpq(3,0)=kpq(3,0)+1
                call fill_indices(col(:),2,1,a,b,i,i,2) 
                kpq(:,cntf)=col(:)
             end do
          end if
       end do
    end do
       
!a=b,i|=j Ne+Mg+  i on Ne, j=n3sMg
  
    do ih=2,ncentre
       i=hcentre(ih)
       do jh=ih+1,ncentre
          j=hcentre(jh)
          isym1=MT(orbSym(i),orbSym(j))
          if (isym1 .eq. nirrep) then
             do ap=nOcc+1,nBas
                a=roccnum(ap)
                cntf=cntf+1
                kpq(4,0)=kpq(4,0)+1
                call fill_indices(col(:),2,1,a,a,i,j,3) 
                kpq(:,cntf)=col(:)
             end do
          end if
       end do
    end do
    
!a|=b,i|=j spin I Ne+Mg+ i on Ne, j=n3sMg

    do ih=2,ncentre
       i=hcentre(ih)
        do jh=ih+1,ncentre
           j=hcentre(jh)
           isym1=MT(orbSym(i),orbSym(j))
           do ap=nOcc+1,nBas
              a=roccnum(ap)
              do bp=ap+1,nBas
                 b=roccnum(bp)
                 isym2=MT(orbSym(a),orbSym(b))
                 if(MT(isym1,isym2) .eq. nirrep) then 
                    cntf=cntf+1
                    kpq(5,0)=kpq(5,0)+1
                    call fill_indices(col(:),2,11,a,b,i,j,4)
                    kpq(:,cntf)=col(:)
                 end if
              end do
           end do
        end do
     end do
    
    kpq(:,cntf+1:cntf+kpq(5,0))=kpq(:,cntf+1-kpq(5,0):cntf)
    kpq(2,cntf+1:cntf+kpq(5,0))=12 
    kpq(7,cntf+1:cnti+kpq(5,0))=5
    
    write(6,100) "Number of 2h-2p |abij> configs in the FS ADC", cntf+kpq(5,0)-kpq(1,0)
    write(6,103)
    write(6,*) " Doubly excited |abij> configs allowed in the Fin. St. Manif."
    write(6,101) "CNF","SPN","HL1","HL2","PT1","PT2"
    
    do k=kpq(1,0)+1, cntf+kpq(5,0)
       write(6,102) k,kpq(2,k),kpq(3,k),&
            kpq(4,k),kpq(5,k),kpq(6,k)
    end do
    
    
  end subroutine select_hetdia_au_d

!!$-------------------------------------------------------------

  subroutine get_ncnfi_ryd(kpq,ncnfi,nryd)

!!$ Subroutine finds the maximum number of Rydberg states in a given basis. 
!!$ It is used to identify the Rydberg series of Fano states (initial states).
    
    integer, dimension(nBas), intent(out) :: ncnfi
    integer, dimension(7,0:nBas**2*4*nOcc**2),intent(in) :: kpq
    integer, intent(out) :: nryd
    
    integer :: i,j,a,cnt

    cnt=0
    
    do i=1,kpq(1,0)
       j=kpq(3,i)
       if((j .eq. hinit)) then 
          cnt=cnt+1
          ncnfi(cnt)=i
       end if
    end do
    nryd=cnt
    
    write(6,*) 'Number of rydberg states', nryd
    write(6,*) 'Cnfs contributing to the Rydbergs', ncnfi(1:nryd)
    
  end subroutine get_ncnfi_ryd

!!$---------------------------------------------
!!$---------------------------------------------

  subroutine select_fstate_ryd(ndim,stsel,vec,nstate,nisri)

    integer, intent(in) :: ndim
    real(d), intent(in) :: stsel
    real(d), dimension(ndim),intent(in) :: vec 
    integer, intent(out) :: nstate
    integer, dimension(ndim),intent(out) :: nisri

    integer :: i
    integer, dimension(ndim) :: indx
    
    nstate=0
    
!!$    call dsortqx("D",ndim,vec(:),1,indx(:))
    
    call dsortindxa1("D",ndim,vec(:),indx(:))
    
    nisri(:)=-1
    
    do i=1,ndim
    if(vec(indx(i)) .ge. stsel) then
       nstate=nstate+1
       nisri(i)=indx(i)
       write(6,*) "selecting states nr.",indx(i)
    else
       exit
    end if
    end do
    
  end subroutine select_fstate_ryd

!!$---------------------------------------------------
!!$---------------------------------------------------

  subroutine select_fstate(ndim,stsel,vec,nstate,nisri)

    integer, intent(in) :: ndim
    real(d), intent(in) :: stsel
    real(d), dimension(ndim),intent(in) :: vec 
    integer,intent(out) :: nstate
    integer, dimension(ndim),intent(out) :: nisri

    integer :: i
    integer, dimension(ndim) :: indx
    real(d), dimension(ndim) :: coeff
    
    nstate=0
    
    do i=1,ndim
       coeff(i)=abs(vec(i))**2 
    end do

!!$    call dsortqx("D",ndim,coeff(:),1,indx(:))
   
    call dsortindxa1("D",ndim,coeff(:),indx(:)) 
    
    nisri(:)=-1
    
    do i=1,ndim
    if(coeff(indx(i)) .ge. stsel) then
       nstate=nstate+1
       nisri(i)=indx(i)
    else
       exit
    end if
    end do
    
    
  end subroutine select_fstate
    
!!$-----------------------------------------------

end module select_fano

    
    
