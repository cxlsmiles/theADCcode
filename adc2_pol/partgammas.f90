module partgammas
  
  use constants
  use parameters
  use get_matrix
  use adc_ph
  use misc
  
  implicit none
  
contains
!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------
!----------------ADC2E ROUTINES------------------------------------------------
!------------------------------------------------------------------------------

   subroutine getf_HIJ_adc2e(ndimi,ndimf,nstate,nstate1,kpq_init,kpq_fin,&
              fstate,fen,arr_fin,gammavec,hpdimf,newnstate,fsen,nfsen,ninista,tvecs)

! fstate and fen enthaölten die eigenvektoren und energien des initial state
     
     integer, intent(in) :: ndimi,ndimf,nstate,nstate1,hpdimf,ninista
     integer, intent(out) :: newnstate 
     integer, dimension(7,0:nBas**2*4*nOcc**2),intent(in) :: kpq_init,kpq_fin
     
     real(d), dimension(ndimf,nstate1),intent(in):: arr_fin
!     real(d), dimension(ndimi),intent(in):: fstate
     real(d), dimension(ndimi,nstate),intent(in):: fstate
     real(d), dimension(nstate),intent(in) :: fen
     real(d), dimension(nstate1),intent(in) :: fsen
     real(d), dimension(nstate1),intent(out) :: nfsen
     real(d), dimension(ndimf,1000), intent(inout) :: tvecs
     
     real(d) ::  matel, pgamma,ea,eb,ej,ek,cntr,etest
     integer :: a,b,j,k,a1,b1,j1,k1,spin,spin1,type,type1,pos,meth
     real(d), dimension(nstate1,nstate),intent(out) :: gammavec
     
     integer :: ic,jc,kc,lc,lim1i,lim2i,lim1j,lim2j,cnt,i,lccnt
     real(d), dimension(ndimf) :: tempvec1,tempvec2
     logical :: diag

     meth=21!method 21-adc2e

!     do kc=1,nstate !counts available initial states
        write(6,*) 'Doing gammas for state nr.',ninista,fen(ninista)
        
        cnt=0
        pgamma=0._d
           
        tempvec1(:)=0._d
        tempvec2(:)=0._d

        
        do jc=1,ndimf !counts final cnfs
           call get_indices1(kpq_fin(:,jc),spin1,a1,b1,j1,k1,type1)
!            write(6,*) 'Outer Final jc',jc
           
           do ic=1,ndimi !counts initial cnfs
              call get_indices1(kpq_init(:,ic),spin,a,b,j,k,type)
!             write(6,*) 'Initial ic',ic
!              write(6,*) 't,t1',type,type1
!              write(6,*) 'a,j',a,j
!              write(6,*) 'a1,j1',a1,j1
            
             
              
              pos=6*type+type1
              
!!$ For historical reasons the confs are ordered in the following blocks: 
!!$ 0-single, 1-aajj,2-abjj,3-aajk,4-abjk(1),5-abjk(2)
              
              select case (pos)
              case(0)!0,0 block
                 call adc2ext_0_0(meth,a,j,a1,j1,matel)
              case(1)!0,1 block
                 matel=C5_ph_2p2h(a,j,a1,j1)
              case(2)!0,2 block
                 matel=C4_ph_2p2h(a,j,a1,b1,j1)
              case(3)!0,3 block
                 matel=C3_ph_2p2h(a,j,a1,j1,k1)
              case(4)!0,4 block
                 matel=C1_ph_2p2h(a,j,a1,b1,j1,k1)
              case(5)!0,5 block
                 matel=C2_ph_2p2h(a,j,a1,b1,j1,k1)
              case(6)!1,0 block
                 matel=C5_ph_2p2h(a1,j1,a,j)
              case(7)!1,1 block
                 ea=e(a)
                 ej=e(j)
                 call adc2ext_1_1(meth,a,j,a1,j1,ea,ej,matel)
              case(8)!1,2 block
                 matel=C_2_1(a1,b1,j1,a,j)
              case(9)!1,3 block
                 matel=C_3_1(a1,j1,k1,a,j)
              case(10)!1,4 block
                 matel=C_4i_1(a1,b1,j1,k1,a,j)
              case(11)!1,5 block
                 matel=C_4ii_1(a1,b1,j1,k1,a,j)
              case(12)!2,0 block
                 matel=C4_ph_2p2h(a1,j1,a,b,j)
              case(13)!2,1 block
                 matel=C_2_1(a,b,j,a1,j1)
              case(14)!2,2 block
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 call adc2ext_2_2(meth,a,b,j,a1,b1,j1,ea,eb,ej,matel)
              case(15)!2,3 block
                 matel=C_3_2(a1,j1,k1,a,b,j)
              case(16)!2,4 block
                 matel=C_4i_2(a1,b1,j1,k1,a,b,j)
              case(17)!2,5
                 matel=C_4ii_2(a1,b1,j1,k1,a,b,j)
              case(18)!3,0
                 matel=C3_ph_2p2h(a1,j1,a,j,k)
              case(19)!3,1
                 matel=C_3_1(a,j,k,a1,j1)
              case(20)!3,2
                 matel=C_3_2(a,j,k,a1,b1,j1)
              case(21)!3,3
                 ea=e(a)
                 ek=e(k)
                 ej=e(j)
                 call adc2ext_3_3(meth,a,j,k,a1,j1,k1,ea,ej,ek,matel)
              case(22)!3,4
                 matel=C_4i_3(a1,b1,j1,k1,a,j,k)
              case(23)!3,5
                 matel=C_4ii_3(a1,b1,j1,k1,a,j,k)
              case(24)!4,0
                 matel=C1_ph_2p2h(a1,j1,a,b,j,k) 
              case(25)
                 matel=C_4i_1(a,b,j,k,a1,j1)
              case(26)!4,2
                 matel=C_4i_2(a,b,j,k,a1,b1,j1)
              case(27)!4,3
                 matel=C_4i_3(a,b,j,k,a1,j1,k1)
              case(28)!4,4
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call  adc2ext_4i_4i(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)
              case(29)!4,5
                 matel=C_4i_4ii(a,b,j,k,a1,b1,j1,k1)
              case(30)!5,0
                 matel=C2_ph_2p2h(a1,j1,a,b,j,k)
              case(31)!5,1
                 matel=C_4ii_1(a,b,j,k,a1,j1)
              case(32)!5,2
                 matel=C_4ii_2(a,b,j,k,a1,b1,j1)
              case(33)!5,3
                 matel=C_4ii_3(a,b,j,k,a1,j1,k1)
              case(34)!5,4
                 matel=C_4i_4ii(a1,b1,j1,k1,a,b,j,k)
              case(35)!5,5
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call adc2ext_4ii_4ii(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)   
              end select

              diag=(a .eq. a1).and.(b .eq. b1).and.(j .eq. j1).and.(k .eq. k1).and.(type .eq. type1)
              if(diag) then
                 cnt=cnt+1
                 matel=matel-fen(ninista)
              end if

!              write(6,*) 'Ninista Energie Test',fen(ninista)

              tempvec1(jc)=tempvec1(jc)+matel*fstate(ic,ninista) !(ninista,nstate)
           end do
        end do

!         do lc=1,100
!         write(6,*) 'tv1',tempvec1(lc)
!         end do
        
        write(6,*) 'hpdimf',hpdimf
        write(6,*) 'back nstate1',nstate1
        lccnt=0
        cntr=0.0
        do lc=1,nstate1 !counts available final states (nstate viele in Main)
           tempvec2(:)=arr_fin(:,lc)

             do i=1,hpdimf
               cntr=cntr+tempvec2(i)**2
             end do
!             write(6,*) 'cntr',cntr
              etest=fsen(lc)
!            if(etest.lt.fen(ninista)) then
!             if(etest.lt.0.5) then
!              lccnt=lccnt+1
!              tvecs(:,lccnt)=tempvec2(:)
!              pgamma=dsp(ndimf,tempvec1(:),tempvec2(:))
!              gammavec(lccnt,1)=2._d*pi*pgamma**2        ! geht nur für lc1cnt << nlim1
!              nfsen(lccnt)=fsen(lc)
!            else
             if(cntr.gt.0.2) then  ! Parameter to select the weigt of the contributing final eigencevtors (states)
             lccnt=lccnt+1
             tvecs(:,lccnt)=tempvec2(:)
             pgamma=dsp(ndimf,tempvec1(:),tempvec2(:))
             gammavec(lccnt,1)=2._d*pi*pgamma**2
             nfsen(lccnt)=fsen(lc)
             end if
!            end if
             cntr=0.0
        end do
!     end do
       write(6,*) 'lccnt',lccnt
       newnstate=lccnt
     
     write(6,*) "Warning. Non-orthogonal iniitial and final subspaces", cnt

   end subroutine getf_HIJ_adc2e

!!$-------------------------------------------------------------------------------------
   subroutine read_HIJ_adc2e(ndimi,ndimf,nstatei,nstatef,kpq_init,kpq_fin,fanostate,fen,gammavec)
     
     integer, intent(in) :: ndimi,ndimf,nstatei,nstatef 
     integer, dimension(7,0:nBas**2*4*nOcc**2),intent(in) :: kpq_init,kpq_fin
     
     real(d), dimension(ndimi,nstatei),intent(in):: fanostate
     real(d), dimension(nstatei),intent(in) :: fen
     
     real(d) ::  matel,pgamma,ea,eb,ej,ek
     integer :: a,b,j,k,a1,b1,j1,k1,spin,spin1,type,type1,pos,meth
     real(d), dimension(nstatef,nstatei),intent(out) :: gammavec
     
     real(d), dimension(:,:), allocatable :: rvec
     integer :: ic,jc,kc,lc1,lc2,lim1i,lim2i,lim1j,lim2j,cnt,mmr,nlim1,nlim2,fnm,nvecout
     real(d), dimension(ndimf) :: tempvec1,tempvec2
     logical :: diag

     mmr=32768000!corr. to 250 Mb mem used in reading Lan. vecs from the saved file
     meth=21!method 21-adc2e
     
     nlim1=mmr/ndimf
     nlim2=lancstates/nlim1+1

     allocate(rvec(ndimf,nlim1))

     do kc=1,nstatei !counts available initial states
        write(6,*) 'Doing gammas for the initial state nr.',kc
        
        cnt=0
        pgamma=0._d
           
        tempvec1(:)=0._d
        tempvec2(:)=0._d

        
        do jc=1,ndimf !counts final cnfs
           call get_indices1(kpq_fin(:,jc),spin1,a1,b1,j1,k1,type1)
           
           do ic=1,ndimi !counts initial cnfs
              call get_indices1(kpq_init(:,ic),spin,a,b,j,k,type)
            
              pos=6*type+type1
              
!!$ For historical reasons the confs are ordered in the following blocks: 
!!$ 0-single, 1-aajj,2-abjj,3-aajk,4-abjk(1),5-abjk(2)
              
              select case (pos)
              case(0)!0,0 block
                 call adc2ext_0_0(meth,a,j,a1,j1,matel)
              case(1)!0,1 block
                 matel=C5_ph_2p2h(a,j,a1,j1)
              case(2)!0,2 block
                 matel=C4_ph_2p2h(a,j,a1,b1,j1)
              case(3)!0,3 block
                 matel=C3_ph_2p2h(a,j,a1,j1,k1)
              case(4)!0,4 block
                 matel=C1_ph_2p2h(a,j,a1,b1,j1,k1)
              case(5)!0,5 block
                 matel=C2_ph_2p2h(a,j,a1,b1,j1,k1)
              case(6)!1,0 block
                 matel=C5_ph_2p2h(a1,j1,a,j)
              case(7)!1,1 block
                 ea=e(a)
                 ej=e(j)
                 call adc2ext_1_1(meth,a,j,a1,j1,ea,ej,matel)
              case(8)!1,2 block
                 matel=C_2_1(a1,b1,j1,a,j)
              case(9)!1,3 block
                 matel=C_3_1(a1,j1,k1,a,j)
              case(10)!1,4 block
                 matel=C_4i_1(a1,b1,j1,k1,a,j)
              case(11)!1,5 block
                 matel=C_4ii_1(a1,b1,j1,k1,a,j)
              case(12)!2,0 block
                 matel=C4_ph_2p2h(a1,j1,a,b,j)
              case(13)!2,1 block
                 matel=C_2_1(a,b,j,a1,j1)
              case(14)!2,2 block
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 call adc2ext_2_2(meth,a,b,j,a1,b1,j1,ea,eb,ej,matel)
              case(15)!2,3 block
                 matel=C_3_2(a1,j1,k1,a,b,j)
              case(16)!2,4 block
                 matel=C_4i_2(a1,b1,j1,k1,a,b,j)
              case(17)!2,5
                 matel=C_4ii_2(a1,b1,j1,k1,a,b,j)
              case(18)!3,0
                 matel=C3_ph_2p2h(a1,j1,a,j,k)
              case(19)!3,1
                 matel=C_3_1(a,j,k,a1,j1)
              case(20)!3,2
                 matel=C_3_2(a,j,k,a1,b1,j1)
              case(21)!3,3
                 ea=e(a)
                 ek=e(k)
                 ej=e(j)
                 call adc2ext_3_3(meth,a,j,k,a1,j1,k1,ea,ej,ek,matel)
              case(22)!3,4
                 matel=C_4i_3(a1,b1,j1,k1,a,j,k)
              case(23)!3,5
                 matel=C_4ii_3(a1,b1,j1,k1,a,j,k)
              case(24)!4,0
                 matel=C1_ph_2p2h(a1,j1,a,b,j,k) 
              case(25)
                 matel=C_4i_1(a,b,j,k,a1,j1)
              case(26)!4,2
                 matel=C_4i_2(a,b,j,k,a1,b1,j1)
              case(27)!4,3
                 matel=C_4i_3(a,b,j,k,a1,j1,k1)
              case(28)!4,4
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call  adc2ext_4i_4i(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)
              case(29)!4,5
                 matel=C_4i_4ii(a,b,j,k,a1,b1,j1,k1)
              case(30)!5,0
                 matel=C2_ph_2p2h(a1,j1,a,b,j,k)
              case(31)!5,1
                 matel=C_4ii_1(a,b,j,k,a1,j1)
              case(32)!5,2
                 matel=C_4ii_2(a,b,j,k,a1,b1,j1)
              case(33)!5,3
                 matel=C_4ii_3(a,b,j,k,a1,j1,k1)
              case(34)!5,4
                 matel=C_4i_4ii(a1,b1,j1,k1,a,b,j,k)
              case(35)!5,5
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call adc2ext_4ii_4ii(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)   
              end select

              diag=(a .eq. a1).and.(b .eq. b1).and.(j .eq. j1).and.(k .eq. k1).and.(type .eq. type1)
              if(diag) then
                 cnt=cnt+1
                 matel=matel-fen(kc)
              end if
              
              tempvec1(jc)=tempvec1(jc)+matel*fanostate(ic,kc)
           end do
        end do

        do lc2=1,nlim2 !counts the batches of read final states
           fnm=2
           call readvct1(ndimf,fnm,(lc2-1)*nlim1+1,lc2*nlim1,rvec,nvecout)
           do lc1=1,nvecout !counts vectors inside each batch
              tempvec2(:)=rvec(:,lc1)
              pgamma=dsp(ndimf,tempvec1(:),tempvec2(:))
              gammavec((lc2-1)*nlim1+lc1,kc)=2._d*pi*pgamma**2
           end do
        end do
        write(6,*) ' vectors read', (lc2-2)*nlim1+nvecout

     end do
     
     deallocate(rvec)
     write(6,*) "Warning. Non-orthogonal iniitial and final subspaces", cnt

   end subroutine read_HIJ_adc2e
!!$----------------------------------------------------------------------------
!!$-------------------------------------------------------------------------------------

   subroutine get_mgamma_adc2e(ndimi,ndimf,nstate,kpq_init,kpq_fin,fanostate,fen,mgammavec)
     
     integer, intent(in) :: ndimi,ndimf,nstate
     integer, dimension(7,0:nBas**2*4*nOcc**2),intent(in) :: kpq_init,kpq_fin
     
     real(d), dimension(ndimi,nstate),intent(in):: fanostate
     real(d), dimension(nstate),intent(in) :: fen
     
     real(d) ::  matel,pgamma,ea,eb,ej,ek
     integer :: a,b,j,k,a1,b1,j1,k1,spin,spin1,type,type1,pos,meth
     real(d), dimension(ndimf),intent(out) :: mgammavec
     
     integer :: ic,jc,lim1i,lim2i,lim1j,lim2j,cnt,nvecout
     real(d), dimension(ndimf) :: tempvec1
     logical :: diag

     meth=21!method 21-adc2e
     
     write(6,*) 'Doing gammas for the initial state nr.', ninista
        
        cnt=0
        pgamma=0._d
           
        tempvec1(:)=0._d
        
        do jc=1,ndimf !counts final cnfs
           call get_indices1(kpq_fin(:,jc),spin1,a1,b1,j1,k1,type1)
           
           do ic=1,ndimi !counts initial cnfs
              call get_indices1(kpq_init(:,ic),spin,a,b,j,k,type)
            
              pos=6*type+type1
              
!!$ For historical reasons the confs are ordered in the following blocks: 
!!$ 0-single, 1-aajj,2-abjj,3-aajk,4-abjk(1),5-abjk(2)
              
              select case (pos)
              case(0)!0,0 block
                 call adc2ext_0_0(meth,a,j,a1,j1,matel)
              case(1)!0,1 block
                 matel=C5_ph_2p2h(a,j,a1,j1)
              case(2)!0,2 block
                 matel=C4_ph_2p2h(a,j,a1,b1,j1)
              case(3)!0,3 block
                 matel=C3_ph_2p2h(a,j,a1,j1,k1)
              case(4)!0,4 block
                 matel=C1_ph_2p2h(a,j,a1,b1,j1,k1)
              case(5)!0,5 block
                 matel=C2_ph_2p2h(a,j,a1,b1,j1,k1)
              case(6)!1,0 block
                 matel=C5_ph_2p2h(a1,j1,a,j)
              case(7)!1,1 block
                 ea=e(a)
                 ej=e(j)
                 call adc2ext_1_1(meth,a,j,a1,j1,ea,ej,matel)
              case(8)!1,2 block
                 matel=C_2_1(a1,b1,j1,a,j)
              case(9)!1,3 block
                 matel=C_3_1(a1,j1,k1,a,j)
              case(10)!1,4 block
                 matel=C_4i_1(a1,b1,j1,k1,a,j)
              case(11)!1,5 block
                 matel=C_4ii_1(a1,b1,j1,k1,a,j)
              case(12)!2,0 block
                 matel=C4_ph_2p2h(a1,j1,a,b,j)
              case(13)!2,1 block
                 matel=C_2_1(a,b,j,a1,j1)
              case(14)!2,2 block
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 call adc2ext_2_2(meth,a,b,j,a1,b1,j1,ea,eb,ej,matel)
              case(15)!2,3 block
                 matel=C_3_2(a1,j1,k1,a,b,j)
              case(16)!2,4 block
                 matel=C_4i_2(a1,b1,j1,k1,a,b,j)
              case(17)!2,5
                 matel=C_4ii_2(a1,b1,j1,k1,a,b,j)
              case(18)!3,0
                 matel=C3_ph_2p2h(a1,j1,a,j,k)
              case(19)!3,1
                 matel=C_3_1(a,j,k,a1,j1)
              case(20)!3,2
                 matel=C_3_2(a,j,k,a1,b1,j1)
              case(21)!3,3
                 ea=e(a)
                 ej=e(j)
                 ek=e(k)
                 call adc2ext_3_3(meth,a,j,k,a1,j1,k1,ea,ej,ek,matel)
              case(22)!3,4
                 matel=C_4i_3(a1,b1,j1,k1,a,j,k)
              case(23)!3,5
                 matel=C_4ii_3(a1,b1,j1,k1,a,j,k)
              case(24)!4,0
                 matel=C1_ph_2p2h(a1,j1,a,b,j,k) 
              case(25)
                 matel=C_4i_1(a,b,j,k,a1,j1)
              case(26)!4,2
                 matel=C_4i_2(a,b,j,k,a1,b1,j1)
              case(27)!4,3
                 matel=C_4i_3(a,b,j,k,a1,j1,k1)
              case(28)!4,4
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call  adc2ext_4i_4i(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)
              case(29)!4,5
                 matel=C_4i_4ii(a,b,j,k,a1,b1,j1,k1)
              case(30)!5,0
                 matel=C2_ph_2p2h(a1,j1,a,b,j,k)
              case(31)!5,1
                 matel=C_4ii_1(a,b,j,k,a1,j1)
              case(32)!5,2
                 matel=C_4ii_2(a,b,j,k,a1,b1,j1)
              case(33)!5,3
                 matel=C_4ii_3(a,b,j,k,a1,j1,k1)
              case(34)!5,4
                 matel=C_4i_4ii(a1,b1,j1,k1,a,b,j,k)
              case(35)!5,5
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call adc2ext_4ii_4ii(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)   
              end select

              diag=(a .eq. a1).and.(b .eq. b1).and.(j .eq. j1).and.(k .eq. k1).and.(type .eq. type1)
              if(diag) then
                 cnt=cnt+1
                 matel=matel-fen(ninista)
              end if
              
              tempvec1(jc)=tempvec1(jc)+matel*fanostate(ic,ninista)
           end do
        end do

!!$        do lc2=1,nlim2 !counts the batches of read final states
!!$           fnm=2
!!$           call readvct1(ndimf,fnm,(lc2-1)*nlim1+1,lc2*nlim1,rvec,nvecout)
!!$           do lc1=1,nvecout !counts vectors inside each batch
!!$              tempvec2(:)=rvec(:,lc1)
!!$              pgamma=dsp(ndimf,tempvec1(:),tempvec2(:))
!!$              gammavec((lc2-1)*nlim1+lc1,kc)=2._d*pi*pgamma**2
!!$           end do
!!$        end do
!!$        write(6,*) ' vectors read', (lc2-2)*nlim1+nvecout

        mgammavec(:)=tempvec1(:)
     
     write(6,*) "Warning. Non-orthogonal iniitial and final subspaces", cnt

   end subroutine get_mgamma_adc2e

!!$-------------------------------------------------------------------------------------
!!$---------------------------ADC2 ROUTINES---------------------------------------------
!!$-------------------------------------------------------------------------------------

   subroutine getf_HIJ_adc2(ndimi,ndimf,istate,nistate,kpq_init,kpq_fin,fen,arr_ini,gammavec)
     
     integer, intent(in) :: ndimi,ndimf,istate,nistate
     integer, dimension(7,0:nBas**2*4*nOcc**2),intent(in) :: kpq_init,kpq_fin
     
     real(d), dimension(ndimi,nistate),intent(in):: arr_ini
     real(d), dimension(nistate),intent(in) :: fen
     
     real(d) ::  matel,pgamma,ea,eb,ej,ek
     integer :: a,b,j,k,a1,b1,j1,k1,spin,spin1,type,type1,pos,meth
     real(d), dimension(ndimf),intent(out) :: gammavec
     
     
     integer :: ic,jc,kc,lc,lim1i,lim2i,lim1j,lim2j,cnt
     real(d), dimension(ndimf) :: tempvec1
     logical :: diag

     meth=2!method 2-adc

        write(6,*) 'Doing gammas for state nr.',istate
        
        cnt=0
        pgamma=0._d
           
        tempvec1(:)=0._d

        do jc=1,ndimf !counts final cnfs
           call get_indices1(kpq_fin(:,jc),spin1,a1,b1,j1,k1,type1)
           
           do ic=1,ndimi !counts initial cnfs
              call get_indices1(kpq_init(:,ic),spin,a,b,j,k,type)
            
             
              
              pos=6*type+type1
              
!!$ For historical reasons the confs are ordered in the following blocks: 
!!$ 0-single, 1-aajj,2-abjj,3-aajk,4-abjk(1),5-abjk(2)
              
              matel=0._d
              
              select case (pos)
              case(0)!0,0 block
                 call adc2ext_0_0(meth,a,j,a1,j1,matel)
              case(1)!0,1 block
                 matel=C5_ph_2p2h(a,j,a1,j1)
              case(2)!0,2 block
                 matel=C4_ph_2p2h(a,j,a1,b1,j1)
              case(3)!0,3 block
                 matel=C3_ph_2p2h(a,j,a1,j1,k1)
              case(4)!0,4 block
                 matel=C1_ph_2p2h(a,j,a1,b1,j1,k1)
              case(5)!0,5 block
                 matel=C2_ph_2p2h(a,j,a1,b1,j1,k1)
              case(6)!1,0 block
                 matel=C5_ph_2p2h(a1,j1,a,j)
              case(7)!1,1 block
                 ea=e(a)
                 ej=e(j)
                 call adc2ext_1_1(meth,a,j,a1,j1,ea,ej,matel)
              case(12)!2,0 block
                 matel=C4_ph_2p2h(a1,j1,a,b,j)
              case(14)!2,2 block
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 call adc2ext_2_2(meth,a,b,j,a1,b1,j1,ea,eb,ej,matel)
              case(18)!3,0
                 matel=C3_ph_2p2h(a1,j1,a,j,k)
              case(21)!3,3
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 call adc2ext_3_3(meth,a,j,k,a1,j1,k1,ea,ej,ek,matel)
              case(24)!4,0
                 matel=C1_ph_2p2h(a1,j1,a,b,j,k) 
              case(28)!4,4
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call  adc2ext_4i_4i(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)
              case(30)!5,0
                 matel=C2_ph_2p2h(a1,j1,a,b,j,k)
              case(35)!5,5
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call adc2ext_4ii_4ii(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)   
              end select

              diag=(a .eq. a1).and.(b .eq. b1).and.(j .eq. j1).and.(k .eq. k1).and.(type .eq. type1)
              if(diag) then
                 cnt=cnt+1
                 matel=matel-fen(istate)
              end if
              
              tempvec1(jc)=tempvec1(jc)+matel*arr_ini(ic,istate)
           end do
        end do

        gammavec(:)=tempvec1(:)
     
     write(6,*) "Warning. Non-orthogonal iniitial and final subspaces", cnt

   end subroutine getf_HIJ_adc2
!!$------------------------------------------------------
!!$-------------------------------------------------------------------------------------
!!$------------------------------------------------------
   subroutine read_HIJ_adc2(ndimi,ndimf,nstate,nstate1,kpq_init,kpq_fin,fstate,fen,gammavec)
     
     integer, intent(in) :: ndimi,ndimf,nstate,nstate1  
     integer, dimension(7,0:nBas**2*4*nOcc**2),intent(in) :: kpq_init,kpq_fin
     
     real(d), dimension(ndimi,nstate),intent(in):: fstate
     real(d), dimension(nstate),intent(in) :: fen
     
     real(d) ::  matel,pgamma,ea,eb,ej,ek
     integer :: a,b,j,k,a1,b1,j1,k1,spin,spin1,type,type1,pos,meth
     real(d), dimension(nstate1,nstate),intent(out) :: gammavec
     
     real(d), dimension(:,:), allocatable :: rvec 
     integer :: ic,jc,kc,lc1,lc2,lim1i,lim2i,lim1j,lim2j,cnt,mmr,nlim1,nlim2,fnm,nvecout
     real(d), dimension(ndimf) :: tempvec1,tempvec2
     logical :: diag

     mmr=32768000!corr. to 250 Mb mem used in reading Lan. vecs from the saved file
     meth=2!method 2-adc

     nlim1=mmr/ndimf
     nlim2=lancstates/nlim1+1
     allocate(rvec(ndimf,nlim1))

     do kc=1,nstate !counts available initial states
        write(6,*) 'Doing gammas for state nr.',kc
        
        cnt=0
        pgamma=0._d
           
        tempvec1(:)=0._d
        tempvec2(:)=0._d

        
        do jc=1,ndimf !counts final cnfs
           call get_indices1(kpq_fin(:,jc),spin1,a1,b1,j1,k1,type1)
           
           do ic=1,ndimi !counts initial cnfs
              call get_indices1(kpq_init(:,ic),spin,a,b,j,k,type)
            
             
              
              pos=6*type+type1
              
!!$ For historical reasons the confs are ordered in the following blocks: 
!!$ 0-single, 1-aajj,2-abjj,3-aajk,4-abjk(1),5-abjk(2)
              
              matel=0._d
              
              select case (pos)
              case(0)!0,0 block
                 call adc2ext_0_0(meth,a,j,a1,j1,matel)
              case(1)!0,1 block
                 matel=C5_ph_2p2h(a,j,a1,j1)
              case(2)!0,2 block
                 matel=C4_ph_2p2h(a,j,a1,b1,j1)
              case(3)!0,3 block
                 matel=C3_ph_2p2h(a,j,a1,j1,k1)
              case(4)!0,4 block
                 matel=C1_ph_2p2h(a,j,a1,b1,j1,k1)
              case(5)!0,5 block
                 matel=C2_ph_2p2h(a,j,a1,b1,j1,k1)
              case(6)!1,0 block
                 matel=C5_ph_2p2h(a1,j1,a,j)
              case(7)!1,1 block
                 ea=e(a)
                 ej=e(j)
                 call adc2ext_1_1(meth,a,j,a1,j1,ea,ej,matel)
              case(12)!2,0 block
                 matel=C4_ph_2p2h(a1,j1,a,b,j)
              case(14)!2,2 block
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 call adc2ext_2_2(meth,a,b,j,a1,b1,j1,ea,eb,ej,matel)
              case(18)!3,0
                 matel=C3_ph_2p2h(a1,j1,a,j,k)
              case(21)!3,3
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 call adc2ext_3_3(meth,a,j,k,a1,j1,k1,ea,ej,ek,matel)
              case(24)!4,0
                 matel=C1_ph_2p2h(a1,j1,a,b,j,k) 
              case(28)!4,4
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call  adc2ext_4i_4i(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)
              case(30)!5,0
                 matel=C2_ph_2p2h(a1,j1,a,b,j,k)
              case(35)!5,5
                 ea=e(a)
                 eb=e(b)
                 ej=e(j)
                 ek=e(k)
                 call adc2ext_4ii_4ii(meth,a,b,j,k,a1,b1,j1,k1,ea,eb,ej,ek,matel)   
              end select

              diag=(a .eq. a1).and.(b .eq. b1).and.(j .eq. j1).and.(k .eq. k1).and.(type .eq. type1)
              if(diag) then
                 cnt=cnt+1
                 matel=matel-fen(kc)
              end if
              
              tempvec1(jc)=tempvec1(jc)+matel*fstate(ic,kc)
           end do
        end do

        do lc2=1,nlim2 !counts the batches of read final states
           fnm=2
           call readvct1(ndimf,fnm,(lc2-1)*nlim1+1,lc2*nlim1,rvec,nvecout)
           do lc1=1,nvecout !counts vectors inside each batch
              tempvec2(:)=rvec(:,lc1)
              pgamma=dsp(ndimf,tempvec1(:),tempvec2(:))
              gammavec((lc2-1)*nlim1+lc1,kc)=2._d*pi*pgamma**2
           end do
        end do
        write(6,*) ' vectors read', (lc2-2)*nlim1+nvecout

     end do

     deallocate(rvec)
     write(6,*) "Warning. Non-orthogonal iniitial and final subspaces", cnt

   end subroutine read_HIJ_adc2

!!$-------------------------------------------------------------------------------------

   subroutine getf_HIJ_tda(ndimi,ndimf,nstate,nstate1,kpq_init,kpq_fin,fstate,fen,arr_fin,gammavec)
     
     integer, intent(in) :: ndimi,ndimf,nstate,nstate1  
     integer, dimension(7,0:nBas**2*4*nOcc**2),intent(in) :: kpq_init,kpq_fin

     real(d), dimension(ndimf,nstate1),intent(in):: arr_fin
     real(d), dimension(ndimi,nstate),intent(in):: fstate
     real(d), dimension(nstate),intent(in) :: fen
     
     real(d) :: matel, pgamma,ea,eb,ej,ek
     integer :: a,b,j,k,a1,b1,j1,k1,spin,spin1,type,type1,pos,meth
     real(d), dimension(nstate1,nstate),intent(out) :: gammavec
     
     
     integer :: ic,jc,kc,lc,lim1i,lim2i,lim1j,lim2j,cnt
     real(d), dimension(ndimf) :: tempvec1,tempvec2
     logical :: diag

     meth=1!method 1-tda

     do kc=1,nstate !counts available initial states
        write(6,*) 'Doing gammas for state nr.',kc
        
        cnt=0
        pgamma=0._d
           
        tempvec1(:)=0._d
        tempvec2(:)=0._d

        
        do jc=1,ndimf !counts final cnfs
           call get_indices1(kpq_fin(:,jc),spin1,a1,b1,j1,k1,type1)
           
           do ic=1,ndimi !counts initial cnfs
              call get_indices1(kpq_init(:,ic),spin,a,b,j,k,type)
            
             
              
              pos=6*type+type1
              if (pos .ne. 0) then
                 write(6,*) 'Wrong configuration'
                 stop
              end if
              
!!$ For historical reasons the confs are ordered in the following blocks: 
!!$ 0-single, 1-aajj,2-abjj,3-aajk,4-abjk(1),5-abjk(2)
              

              call adc2ext_0_0(meth,a,j,a1,j1,matel)

              diag=(a .eq. a1).and.(b .eq. b1).and.(j .eq. j1).and.(k .eq. k1).and.(type .eq. type1)
              if(diag) then
                 cnt=cnt+1
                 matel=matel-fen(kc)
              end if
              
              tempvec1(jc)=tempvec1(jc)+matel*fstate(ic,kc)
           end do
        end do

        do lc=1,nstate1 !counts available final states
           
           tempvec2(:)=arr_fin(:,lc)
           pgamma=dsp(ndimf,tempvec1(:),tempvec2(:))
           gammavec(lc,kc)=2._d*pi*pgamma**2
           
        end do
     end do
     
     write(6,*) "Warning. Non-orthogonal iniitial and final subspaces", cnt

   end subroutine getf_HIJ_tda

  
   
end module partgammas
     
