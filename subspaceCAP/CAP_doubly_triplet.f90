!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_isr_print_triplet(eigvec1,eigvec2,mat)
use ADC
use def_config
use CAP_MO
implicit none

 type(ADC_eigenvectors), intent(in) :: eigvec1, eigvec2
 double precision, intent(out) :: mat

 double precision :: mat_spac
 integer :: i, j

 mat = 0d0

 open(unit=101,file='ISRmat.tmp')
 open(unit=102,file='ISRmat.txt')

 do i = 1, eigvec1%nconfig
!   do j = 1, eigvec2%nconfig
   do j = 1, i!, eigvec2%nconfig

    call CAP_doubly_spin_adapted_config_triplet(eigvec1%cfg(i),eigvec2%cfg(j),mat_spac)
      mat = mat + eigvec1%coeff(i)*eigvec2%coeff(j)*mat_spac

!      write(101,'(2i4,4X,6(f20.16))')i-1,j-1,eigvec1%coeff(i),eigvec2%coeff(j),mat_spac
!      write(*,'(2i4,4X,1(f20.16))')i-1,j-1,mat_spac
      write(102,'(2(i2,1X),1X,6(f18.12))')i-1,j-1,mat_spac
      write(101,'((f11.8))')mat_spac

   enddo
!      write(100,*)
 enddo

 close(101)
 close(102)

end subroutine test_isr_print_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ADCCAP_DOUBLY_TRIPLET(eigvec1,eigvec2)
use ADC
use def_config
use CAP_MO
implicit none

 type(ADC_eigenvectors), intent(in) :: eigvec1, eigvec2

 double precision :: mat_spac
 integer :: i, j

 do i = 1, eigvec1%nconfig
   do j = 1, eigvec2%nconfig

    call CAP_doubly_spin_adapted_config_triplet(eigvec1%cfg(i),eigvec2%cfg(j),mat_spac)
    HCAPADC(i,j) = mat_spac

   enddo
 enddo

end subroutine ADCCAP_DOUBLY_TRIPLET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CAP_doubly_triplet(eigvec1,eigvec2,mat)
use ADC
use def_config
use CAP_MO
implicit none

 type(ADC_eigenvectors), intent(in) :: eigvec1, eigvec2
 double precision, intent(out) :: mat

 double precision :: mat_spac
 integer :: i, j

 mat = 0d0

 do i = 1, eigvec1%nconfig
   do j = 1, eigvec2%nconfig
      mat = mat + eigvec1%coeff(i)*eigvec2%coeff(j)*HCAPADC(i,j)
   enddo
 enddo

end subroutine CAP_doubly_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Computes the CAP between two spin adapted config.

subroutine CAP_doubly_spin_adapted_config_triplet(spac1,spac2,mat_spac)
use def_config
implicit none

 type(spin_adapted_config), intent(inout) :: spac1, spac2
 double precision, intent(out) :: mat_spac

 double precision :: CAP_doubly_ij_kl_triplet, CAP_doubly_ij_rklm_triplet, &
 CAP_doubly_rijk_slmn_triplet, CAP_doubly_rjii_slkk_triplet, CAP_doubly_smll_rijk_triplet
 integer :: i, j

    mat_spac = 0d0

! only triplet spin

 if(spac1%nholes==2 .and. spac2%nholes==2) then

    mat_spac = CAP_doubly_ij_kl_triplet(spac1%ihole(1),spac1%ihole(2),spac2%ihole(1),spac2%ihole(2))

 elseif(spac1%nholes==2) then

! nico 16.06.2011
! see Yasen for the minus
    mat_spac = -CAP_doubly_ij_rklm_triplet(spac1%ihole(1),spac1%ihole(2),spac2%ipart(1),&
     spac2%ihole(1),spac2%ihole(2),spac2%ihole(3),spac2%typ)

 elseif(spac2%nholes==2) then

! nico 16.06.2011
! see Yasen for the minus
    mat_spac = -CAP_doubly_ij_rklm_triplet(spac2%ihole(1),spac2%ihole(2),spac1%ipart(1),&
     spac1%ihole(1),spac1%ihole(2),spac1%ihole(3),spac1%typ)

 else

    if(spac1%typ==1 .and. spac2%typ==1) then

        mat_spac = CAP_doubly_rjii_slkk_triplet(spac1%ipart(1),spac1%ihole(1),spac1%ihole(2),&
         spac2%ipart(1),spac2%ihole(1),spac2%ihole(2))

    elseif(spac1%typ==1) then

        mat_spac = CAP_doubly_smll_rijk_triplet(spac1%ipart(1),spac1%ihole(1),spac1%ihole(2),&
           spac2%ipart(1),spac2%ihole(1),spac2%ihole(2),spac2%ihole(3),spac2%typ)

    elseif(spac2%typ==1) then

        mat_spac = CAP_doubly_smll_rijk_triplet(spac2%ipart(1),spac2%ihole(1),spac2%ihole(2),&
           spac1%ipart(1),spac1%ihole(1),spac1%ihole(2),spac1%ihole(3),spac1%typ)

    else

        mat_spac = CAP_doubly_rijk_slmn_triplet(spac1%ipart(1),spac1%ihole(1),spac1%ihole(2),&
           spac1%ihole(3),spac1%typ,spac2%ipart(1),spac2%ihole(1),spac2%ihole(2),spac2%ihole(3),spac2%typ)

    endif

 endif

end subroutine CAP_doubly_spin_adapted_config_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_doubly_ij_kl_triplet(i,j,k,l)
use CAP_MO
use dmat_triplet
implicit none

 double precision :: CAP_doubly_ij_kl_triplet
 integer, intent(in) :: i, j, k, l

  CAP_doubly_ij_kl_triplet = 0d0

 if(i==k .and. j==l) CAP_doubly_ij_kl_triplet = CAP_doubly_ij_kl_triplet + D02
 if(i==k) CAP_doubly_ij_kl_triplet = CAP_doubly_ij_kl_triplet - mat_CAP_MO(l,j)
 if(j==l) CAP_doubly_ij_kl_triplet = CAP_doubly_ij_kl_triplet - mat_CAP_MO(k,i)

 if(i==l) CAP_doubly_ij_kl_triplet = CAP_doubly_ij_kl_triplet + mat_CAP_MO(k,j)
 if(j==k) CAP_doubly_ij_kl_triplet = CAP_doubly_ij_kl_triplet + mat_CAP_MO(l,i)

 CAP_doubly_ij_kl_triplet = CAP_doubly_ij_kl_triplet + 1d0*tmat(i,j,k,l) &
 + d21(i,j,k,l) + 1d0*d22(i,j,k,l) + 1d0*d23(i,j,k,l) &
 + 1d0*d24(i,j,k,l) + 1d0*d25(i,j,k,l) + 1d0*d26(i,j,k,l) &
 + 1d0*d27(i,j,k,l) + 1d0*d28(i,j,k,l) + 1d0*d29(i,j,k,l) &
 + 1d0*d210(i,j,k,l) + 1d0*d211(i,j,k,l) + 1d0*d212(i,j,k,l) &
 + 1d0*d213(i,j,k,l) + 1d0*d214(i,j,k,l) + 1d0*d215(i,j,k,l) &
!! + h.c.
 + 1d0*tmat(k,l,i,j) &
 + 1d0*d21(k,l,i,j) + 1d0*d22(k,l,i,j) + 1d0*d23(k,l,i,j) &
 + 1d0*d24(k,l,i,j) + 1d0*d25(k,l,i,j) + 1d0*d26(k,l,i,j) &
 + 1d0*d27(k,l,i,j) + 1d0*d28(k,l,i,j) + 1d0*d29(k,l,i,j) &
 + 1d0*d210(k,l,i,j) + 1d0*d211(k,l,i,j) + 1d0*d212(k,l,i,j) &
 + 1d0*d213(k,l,i,j) + 1d0*d214(k,l,i,j) + 1d0*d215(k,l,i,j) 

end function CAP_doubly_ij_kl_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_doubly_ij_rklm_triplet(i,j,r,k,l,m,typ)
use CAP_MO
use dmat_triplet
implicit none

 double precision :: CAP_doubly_ij_rklm_triplet
 integer, intent(in) :: i, j, r, k, l, m, typ

 CAP_doubly_ij_rklm_triplet = 0d0

    if(typ == 1) then

! warning k and l are all mixed up compare to pdf file
      if(i==k .and. j==l) CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet - mat_CAP_MO(l,r) 
      if(i==l .and. j==k) CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + mat_CAP_MO(l,r) 

      CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + t21(i,j,r,k,l,m)

    elseif(typ == 21) then

      if(i==l .and. j==m) CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + mat_CAP_MO(k,r) 
      if(i==k .and. j==m) CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet - mat_CAP_MO(l,r) 

      CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + t22(i,j,r,k,l,m)

    elseif(typ == 22) then

      if(i==k .and. j==l) CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + mat_CAP_MO(m,r) 
      if(i==l .and. j==m) CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + mat_CAP_MO(k,r) 

      CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + t23(i,j,r,k,l,m)

    elseif(typ == 23) then

      if(i==k .and. j==l) CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + mat_CAP_MO(m,r) 
      if(i==k .and. j==m) CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet - mat_CAP_MO(l,r) 

      CAP_doubly_ij_rklm_triplet = CAP_doubly_ij_rklm_triplet + t24(i,j,r,k,l,m)

    endif

end function CAP_doubly_ij_rklm_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_doubly_rjii_slkk_triplet(r,j,i,s,l,k)
use CAP_MO
implicit none

 double precision :: CAP_doubly_rjii_slkk_triplet
 integer, intent(in) :: r, i, j, k, s, l

 CAP_doubly_rjii_slkk_triplet = 0d0

  if(r==s .and. j==l .and. i == k) CAP_doubly_rjii_slkk_triplet = CAP_doubly_rjii_slkk_triplet + D00

  if(j==l .and. i == k) CAP_doubly_rjii_slkk_triplet = CAP_doubly_rjii_slkk_triplet + mat_CAP_MO(r,s)
  if(r == s .and. j==l .and. i == k) CAP_doubly_rjii_slkk_triplet = CAP_doubly_rjii_slkk_triplet - 2d0 * mat_CAP_MO(k,i)
  if(r == s .and. j==k .and. i == l) CAP_doubly_rjii_slkk_triplet = CAP_doubly_rjii_slkk_triplet + 1d0 * mat_CAP_MO(k,i)
  if(r == s .and. i==k) CAP_doubly_rjii_slkk_triplet = CAP_doubly_rjii_slkk_triplet - 1d0 * mat_CAP_MO(l,j)

end function CAP_doubly_rjii_slkk_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_doubly_smll_rijk_triplet(s,m,l,r,i,j,k,typ)
use CAP_MO
implicit none

 double precision :: CAP_doubly_smll_rijk_triplet
 integer, intent(in) :: s, m, l, r, i, j, k, typ

 CAP_doubly_smll_rijk_triplet = 0d0

 if(typ==21) then
 
    if(r==s .and. i==m .and. j == l) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet - mat_CAP_MO(l,k)
    if(r==s .and. i==l .and. j == m) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet + mat_CAP_MO(l,k)
    if(r==s .and. i==m .and. k == l) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet - mat_CAP_MO(l,j)
    if(r==s .and. j==m .and. k == l) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet + mat_CAP_MO(l,i)

 elseif(typ==22) then

    if(r==s .and. i==m .and. j == l) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet + mat_CAP_MO(l,k)
    if(r==s .and. i==m .and. k == l) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet + mat_CAP_MO(l,j)
    if(r==s .and. i==l .and. k == m) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet - mat_CAP_MO(l,j)
    if(r==s .and. j==l .and. k == m) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet - mat_CAP_MO(l,i)

 elseif(typ==23) then

    if(r==s .and. i==l .and. j == m) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet - mat_CAP_MO(l,k)
    if(r==s .and. i==l .and. k == m) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet + mat_CAP_MO(l,j)
    if(r==s .and. j==m .and. k == l) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet - mat_CAP_MO(l,i)
    if(r==s .and. j==l .and. k == m) CAP_doubly_smll_rijk_triplet = CAP_doubly_smll_rijk_triplet + mat_CAP_MO(l,i)

 endif

end function CAP_doubly_smll_rijk_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_doubly_rijk_slmn_triplet(r,i,j,k,typ1,s,l,m,n,typ2)
use CAP_MO
implicit none

 double precision :: CAP_doubly_rijk_slmn_triplet
 integer, intent(in) :: r, i, j, k, s, l, m, n, typ1, typ2
 integer :: rt, it, jt, kt

 CAP_doubly_rijk_slmn_triplet = 0d0

 if(typ1==21 .and. typ2==21) then

   if(r==s .and. i==l .and. j == m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + D00
   if(i==l .and. j==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(r,s)
   if(r==s .and. i==l .and. j==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(n,k)
   if(r==s .and. i==l .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(m,j)
   if(r==s .and. i==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(l,j)
   if(r==s .and. j==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(l,i)
   if(r==s .and. j==l .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(m,i)

 elseif(typ1==22 .and. typ2==21) then

   if(r==s .and. i==l .and. j==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(m,k)
   if(r==s .and. i==m .and. j==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(l,k)
   if(r==s .and. i==l .and. k==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(n,j)

 elseif(typ1==23 .and. typ2==21) then

   if(r==s .and. j==l .and. k==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(n,i)

 elseif(typ1==21 .and. typ2==22) then

   if(r==s .and. i==l .and. j==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(m,k)
   if(r==s .and. i==l .and. k==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(n,j)
   if(r==s .and. j==l .and. k==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(n,i)

 elseif(typ1==22 .and. typ2==22) then

   if(r==s .and. i==l .and. j == m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + D00
   if(i==l .and. j==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(r,s)
   if(r==s .and. i==l .and. j==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(n,k)
   if(r==s .and. i==l .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(m,j)
   if(r==s .and. j==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(l,i)

 elseif(typ1==23 .and. typ2==22) then

   if(r==s .and. i==m .and. j==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(l,k)
   if(r==s .and. i==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(l,j)
   if(r==s .and. j==l .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(m,i)

 elseif(typ1==21 .and. typ2==23) then

   if(r==s .and. i==m .and. j==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(l,k)

 elseif(typ1==22 .and. typ2==23) then

   if(r==s .and. i==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(l,j)
   if(r==s .and. j==l .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(m,i)
   if(r==s .and. j==l .and. k==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(n,i)


 elseif(typ1==23 .and. typ2==23) then

   if(r==s .and. i==l .and. j == m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + D00
   if(i==l .and. j==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(r,s)
   if(r==s .and. i==l .and. j==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(n,k)
   if(r==s .and. i==l .and. j==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(m,k)
   if(r==s .and. i==l .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(m,j)
   if(r==s .and. i==l .and. k==m) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet + mat_CAP_MO(n,j)
   if(r==s .and. j==m .and. k==n) CAP_doubly_rijk_slmn_triplet = CAP_doubly_rijk_slmn_triplet - mat_CAP_MO(l,i)

 else
   write(*,*)'error in CAP_doubly_rijk_slmn_triplet, I stop'
   stop
 endif

end function CAP_doubly_rijk_slmn_triplet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
