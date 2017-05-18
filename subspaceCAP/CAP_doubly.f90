!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine test_isr_print(eigvec1,eigvec2,mat)
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
   do j = 1, eigvec2%nconfig
!   do j = 1, i!, eigvec2%nconfig

    call CAP_doubly_spin_adapted_config(eigvec1%cfg(i),eigvec2%cfg(j),mat_spac)
!      mat = mat + eigvec1%coeff(i)*eigvec2%coeff(j)*mat_spac

!      write(101,'(2i4,4X,6(f20.16))')i-1,j-1,eigvec1%coeff(i),eigvec2%coeff(j),mat_spac
!      write(*,'(2i4,4X,1(f20.16))')i-1,j-1,mat_spac
      write(102,'(2(i2,1X),1X,6(f18.12))')i-1,j-1,mat_spac
      write(101,'((f11.8))')mat_spac

   enddo
!      write(100,*)
 enddo

 close(101)
 close(102)

end subroutine test_isr_print

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ADCCAP_DOUBLY(eigvec1,eigvec2)
use ADC
use def_config
use CAP_MO
implicit none

 type(ADC_eigenvectors), intent(in) :: eigvec1, eigvec2

 double precision :: mat_spac
 integer :: i, j

 do i = 1, eigvec1%nconfig
   do j = 1, eigvec2%nconfig

    call CAP_doubly_spin_adapted_config(eigvec1%cfg(i),eigvec2%cfg(j),mat_spac)
    HCAPADC(i,j) = mat_spac

   enddo
 enddo

end subroutine ADCCAP_DOUBLY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CAP_doubly(eigvec1,eigvec2,mat)
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

end subroutine CAP_doubly

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Computes the CAP between two spin adapted config.

subroutine CAP_doubly_spin_adapted_config(spac1,spac2,mat_spac)
use def_config
implicit none

 type(spin_adapted_config), intent(inout) :: spac1, spac2
 double precision, intent(out) :: mat_spac

 double precision :: CAP_doubly_ij_kl, CAP_doubly_ij_rklm, CAP_doubly_rijk_slmn
 integer :: i, j

    mat_spac = 0d0

! only singlet spin

 if(spac1%nholes==2 .and. spac2%nholes==2) then

    mat_spac = CAP_doubly_ij_kl(spac1%ihole(1),spac1%ihole(2),spac2%ihole(1),spac2%ihole(2))

 elseif(spac1%nholes==2) then

! nico 16.06.2011
! see Yasen for the minus
   mat_spac = -CAP_doubly_ij_rklm(spac1%ihole(1),spac1%ihole(2),spac2%ipart(1),&
&spac2%ihole(1),spac2%ihole(2),spac2%ihole(3),spac2%typ)

 elseif(spac2%nholes==2) then

! nico 16.06.2011
! see Yasen for the minus
    mat_spac = -CAP_doubly_ij_rklm(spac2%ihole(1),spac2%ihole(2),spac1%ipart(1),&
&spac1%ihole(1),spac1%ihole(2),spac1%ihole(3),spac1%typ)

 else

    mat_spac = CAP_doubly_rijk_slmn(spac1%ipart(1),spac1%ihole(1),spac1%ihole(2),spac1%ihole(3),&
             &spac1%typ,spac2%ipart(1),spac2%ihole(1),spac2%ihole(2),spac2%ihole(3),spac2%typ)

 endif

end subroutine CAP_doubly_spin_adapted_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_doubly_ij_kl(i,j,k,l)
use CAP_MO
use dmat
implicit none

 double precision :: CAP_doubly_ij_kl
 integer, intent(in) :: i, j, k, l

  CAP_doubly_ij_kl = 0d0

 if(i==j .and. k==l) then

  if(i==k) then
   CAP_doubly_ij_kl = D02 - 2d0*mat_CAP_MO(i,k) 

   CAP_doubly_ij_kl = CAP_doubly_ij_kl + tmat(1,i,j,k,l)  &
  + d21(1,i,j,k,l) + d22(1,i,j,k,l) + d23(1,i,j,k,l) &
  + d24(1,i,j,k,l) + d25(1,i,j,k,l) + d26(1,i,j,k,l) &
  + d27(1,i,j,k,l) + d28(1,i,j,k,l) + d29(1,i,j,k,l) &
  + d210(1,i,j,k,l) + d211(1,i,j,k,l) + d212(1,i,j,k,l) &
  + d213(1,i,j,k,l) + d214(1,i,j,k,l) + d215(1,i,j,k,l) &
! + h.c.
  + d21(1,k,l,i,j) + d22(1,k,l,i,j) + d23(1,k,l,i,j) &
  + d24(1,k,l,i,j) + d25(1,k,l,i,j) + d26(1,k,l,i,j) &
  + d27(1,k,l,i,j) + d28(1,k,l,i,j) + d29(1,k,l,i,j) &
  + d210(1,k,l,i,j) + d211(1,k,l,i,j) + d212(1,k,l,i,j) &
  + d213(1,k,l,i,j) + d214(1,k,l,i,j) + d215(1,k,l,i,j) 
  
  else

   CAP_doubly_ij_kl = CAP_doubly_ij_kl + tmat(1,i,j,k,l)  &
  + d21(1,i,j,k,l) + d22(1,i,j,k,l) + d23(1,i,j,k,l) &
  + d24(1,i,j,k,l) + d25(1,i,j,k,l) + d26(1,i,j,k,l) &
  + d27(1,i,j,k,l) + d28(1,i,j,k,l) + d29(1,i,j,k,l) &
  + d210(1,i,j,k,l) + d211(1,i,j,k,l) + d212(1,i,j,k,l) &
  + d213(1,i,j,k,l) + d214(1,i,j,k,l) + d215(1,i,j,k,l) &
! + h.c.
  + d21(1,k,l,i,j) + d22(1,k,l,i,j) + d23(1,k,l,i,j) &
  + d24(1,k,l,i,j) + d25(1,k,l,i,j) + d26(1,k,l,i,j) &
  + d27(1,k,l,i,j) + d28(1,k,l,i,j) + d29(1,k,l,i,j) &
  + d210(1,k,l,i,j) + d211(1,k,l,i,j) + d212(1,k,l,i,j) &
  + d213(1,k,l,i,j) + d214(1,k,l,i,j) + d215(1,k,l,i,j)

  endif

 elseif(i==j) then 

  if(i==l) CAP_doubly_ij_kl = CAP_doubly_ij_kl - dsqrt(2d0)*mat_CAP_MO(k,i) 
  if(i==k) CAP_doubly_ij_kl = CAP_doubly_ij_kl - dsqrt(2d0)*mat_CAP_MO(l,i)  

  CAP_doubly_ij_kl = CAP_doubly_ij_kl + tmat(2,k,l,i,j) &
  + d21(2,i,j,k,l) + d22(2,i,j,k,l) + d23(2,i,j,k,l) &
  + d24(2,i,j,k,l) + d25(2,i,j,k,l) + d26(2,i,j,k,l) &
  + d27(2,i,j,k,l) + d28(2,i,j,k,l) + d29(2,i,j,k,l) &
  + d210(2,i,j,k,l) + d211(2,i,j,k,l) + d212(2,i,j,k,l) &
  + d213(2,i,j,k,l) + d214(2,i,j,k,l) + d215(2,i,j,k,l) &
! + h.c.
  + d21(2,k,l,i,j) + d22(2,k,l,i,j) + d23(2,k,l,i,j) &
  + d24(2,k,l,i,j) + d25(2,k,l,i,j) + d26(2,k,l,i,j) &
  + d27(2,k,l,i,j) + d28(2,k,l,i,j) + d29(2,k,l,i,j) &
  + d210(2,k,l,i,j) + d211(2,k,l,i,j) + d212(2,k,l,i,j) &
  + d213(2,k,l,i,j) + d214(2,k,l,i,j) + d215(2,k,l,i,j)

 elseif(k==l) then

  if(k==j) CAP_doubly_ij_kl = CAP_doubly_ij_kl - dsqrt(2d0)*mat_CAP_MO(i,k) 
  if(k==i) CAP_doubly_ij_kl = CAP_doubly_ij_kl - dsqrt(2d0)*mat_CAP_MO(j,k) 

  CAP_doubly_ij_kl = CAP_doubly_ij_kl + tmat(2,i,j,k,l) &
  + d21(2,i,j,k,l) + d22(2,i,j,k,l) + d23(2,i,j,k,l) &
  + d24(2,i,j,k,l) + d25(2,i,j,k,l) + d26(2,i,j,k,l) &
  + d27(2,i,j,k,l) + d28(2,i,j,k,l) + d29(2,i,j,k,l) &
  + d210(2,i,j,k,l) + d211(2,i,j,k,l) + d212(2,i,j,k,l) &
  + d213(2,i,j,k,l) + d214(2,i,j,k,l) + d215(2,i,j,k,l) & 
! + h.c.
  + d21(2,k,l,i,j) + d22(2,k,l,i,j) + d23(2,k,l,i,j) &
  + d24(2,k,l,i,j) + d25(2,k,l,i,j) + d26(2,k,l,i,j) &
  + d27(2,k,l,i,j) + d28(2,k,l,i,j) + d29(2,k,l,i,j) &
  + d210(2,k,l,i,j) + d211(2,k,l,i,j) + d212(2,k,l,i,j) &
  + d213(2,k,l,i,j) + d214(2,k,l,i,j) + d215(2,k,l,i,j)

 else

  if(i==k .and. j==l) CAP_doubly_ij_kl = CAP_doubly_ij_kl + D02   
  if(j==l) CAP_doubly_ij_kl = CAP_doubly_ij_kl - mat_CAP_MO(k,i)  
  if(i==k) CAP_doubly_ij_kl = CAP_doubly_ij_kl - mat_CAP_MO(l,j)  
  if(j==k) CAP_doubly_ij_kl = CAP_doubly_ij_kl - mat_CAP_MO(l,i)
  if(i==l) CAP_doubly_ij_kl = CAP_doubly_ij_kl - mat_CAP_MO(k,j)  

  CAP_doubly_ij_kl = CAP_doubly_ij_kl+ tmat(3,i,j,k,l) 

  CAP_doubly_ij_kl = CAP_doubly_ij_kl &
  + d21(3,i,j,k,l) + d22(3,i,j,k,l) + d23(3,i,j,k,l) &
  + d24(3,i,j,k,l) + d25(3,i,j,k,l) + d26(3,i,j,k,l) &
  + d27(3,i,j,k,l) + d28(3,i,j,k,l) + d29(3,i,j,k,l) &
  + d210(3,i,j,k,l) + d211(3,i,j,k,l) + d212(3,i,j,k,l) &
  + d213(3,i,j,k,l) + d214(3,i,j,k,l) + d215(3,i,j,k,l) &
! + h.c.
  + d21(3,k,l,i,j) + d22(3,k,l,i,j) + d23(3,k,l,i,j) &
  + d24(3,k,l,i,j) + d25(3,k,l,i,j) + d26(3,k,l,i,j) &
  + d27(3,k,l,i,j) + d28(3,k,l,i,j) + d29(3,k,l,i,j) &
  + d210(3,k,l,i,j) + d211(3,k,l,i,j) + d212(3,k,l,i,j) &
  + d213(3,k,l,i,j) + d214(3,k,l,i,j) + d215(3,k,l,i,j)

 endif

end function CAP_doubly_ij_kl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_doubly_ij_rklm(i,j,r,k,l,m,typ)
use CAP_MO
use dmat
implicit none

 double precision :: CAP_doubly_ij_rklm
 integer, intent(in) :: i, j, r, k, l, m, typ

 CAP_doubly_ij_rklm = 0d0

 if(i==j) then

    if(typ == 1) then

      if(i==l) CAP_doubly_ij_rklm = dsqrt(2d0)*mat_CAP_MO(k,r) 

      CAP_doubly_ij_rklm = CAP_doubly_ij_rklm + t21(i,j,r,k,l,m)

    elseif(typ == 21) then

      CAP_doubly_ij_rklm = t23(i,j,r,k,l,m) 

    elseif(typ == 22) then

      CAP_doubly_ij_rklm = t25(i,j,r,k,l,m) 

    endif

 else

    if(typ == 1) then

!old      if(i==k .and. j==l) CAP_doubly_ij_rklm = -mat_CAP_MO(l,r) !&
! new 23.02.2011
      if(i==k .and. j==l) CAP_doubly_ij_rklm = -mat_CAP_MO(l,r)  
      if(i==l .and. j==k) CAP_doubly_ij_rklm = CAP_doubly_ij_rklm - mat_CAP_MO(l,r)  

      CAP_doubly_ij_rklm = CAP_doubly_ij_rklm + t22(i,j,r,k,l,m) 

    elseif(typ == 21) then

      if(j==m .and. i==l) CAP_doubly_ij_rklm = CAP_doubly_ij_rklm - mat_CAP_MO(k,r)/dsqrt(2d0)
      if(j==m .and. i==k) CAP_doubly_ij_rklm = CAP_doubly_ij_rklm + mat_CAP_MO(l,r)*dsqrt(2d0)
      if(j==l .and. i==k) CAP_doubly_ij_rklm = CAP_doubly_ij_rklm - mat_CAP_MO(m,r)/dsqrt(2d0)

      CAP_doubly_ij_rklm = CAP_doubly_ij_rklm + t24(i,j,r,k,l,m) 

    elseif(typ == 22) then

      if(j==m .and. i==l) CAP_doubly_ij_rklm = CAP_doubly_ij_rklm - mat_CAP_MO(k,r)*dsqrt(1.5d0)
      if(j==l .and. i==k) CAP_doubly_ij_rklm = CAP_doubly_ij_rklm + mat_CAP_MO(m,r)*dsqrt(1.5d0)

      CAP_doubly_ij_rklm = CAP_doubly_ij_rklm + t26(i,j,r,k,l,m) 

    endif

 endif

end function CAP_doubly_ij_rklm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_doubly_rijk_slmn(r,i,j,k,typ1,r1,i1,j1,k1,typ2)
use CAP_MO
implicit none

 double precision :: CAP_doubly_rijk_slmn
 integer :: r, i, j, k, r1, i1, j1, k1, typ1, typ2
 integer :: rt, it, jt, kt

! write(*,*)"B",typ1,typ2

 CAP_doubly_rijk_slmn = 0d0

 if(typ1 == 1 .and. typ2 == 1) then

!old
!  if(r==r1 .and. i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + D00
!  if(i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(r,r1)
!  if(r == r1 .and. i==i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 2d0 * mat_CAP_MO(k1,k)
!  if(r == r1 .and. i==j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0 * mat_CAP_MO(i1,j)
!  if(r == r1 .and. j==j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(i1,i)
!  if(r == r1 .and. j==i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0 * mat_CAP_MO(j1,i)

! new 23.02.2011
  if(r==r1 .and. i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + D00
!  WRITE(*,*)"NICO ", D00
  if(i==i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(r,r1)
  if(r == r1 .and. i==i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 2d0 * mat_CAP_MO(j1,j)
  if(r == r1 .and. i==j1 .and. j == i1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0 * mat_CAP_MO(j1,j)
  if(r == r1 .and. j==j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(i1,i)


 elseif(typ1 == 1) then

  if(typ2 == 21) then

!! HERE K1 AND I ARE INTER-CHANGED IN COMPARISON OF YASEN EQUATIONS
!! I SHOULD MAKE SURE THAT IT IS OK

! old
!    if(r==r1 .and. i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - D00
!    if(i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(r,r1)
!    if(r==r1 .and. i == i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 2d0*mat_CAP_MO(j1,j)
!    if(r==r1 .and. i == i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 2d0*mat_CAP_MO(k1,k)
!    if(r==r1 .and. k == k1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(i1,i)
!    if(r==r1 .and. k == i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(k1,i)
!    if(r==r1 .and. i == j1 .and. j == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 4d0*mat_CAP_MO(i1,j)
!    if(r==r1 .and. j == i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 2d0*mat_CAP_MO(j1,i)
!    CAP_doubly_rijk_slmn = 0.5d0*CAP_doubly_rijk_slmn/dsqrt(2d0)

! new 23.02.2011
    if(r1==r .and. i1 == j .and. j1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 2d0*mat_CAP_MO(j,k1)
    if(r1==r .and. i1 == j .and. k1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0*mat_CAP_MO(j,j1)
    if(r1==r .and. j1 == i .and. k1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 2d0*mat_CAP_MO(j,i1)
    if(r1==r .and. j1 == j .and. k1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0*mat_CAP_MO(j,i1)
    if(r1==r .and. i1 == i .and. j1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0*mat_CAP_MO(j,k1)
    if(r1==r .and. i1 == i .and. k1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0*mat_CAP_MO(j,j1)
    CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn/dsqrt(2d0)

  else

!! HERE K1 AND I ARE INTER-CHANGED IN COMPARISON OF YASEN EQUATIONS
!! I SHOULD MAKE SURE THAT IT IS OK

! old
!    if(r==r1 .and. i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - D00
!    if(i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(r,r1)
!    if(r==r1 .and. i == i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 2d0*mat_CAP_MO(j1,j)
!    if(r==r1 .and. i == i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 2d0*mat_CAP_MO(k1,k)
!    if(r==r1 .and. j == i1 .and. k == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(k1,i)
!    if(r==r1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(i1,i)
!    CAP_doubly_rijk_slmn = 0.5d0*dsqrt(1.5d0)*CAP_doubly_rijk_slmn

! new 23.02.2011
    if(r1==r .and. i1 == i .and. j1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(j,k1)
    if(r1==r .and. i1 == i .and. k1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(j,j1)
    if(r1==r .and. i1 == j .and. k1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(j,j1)
    if(r1==r .and. j1 == j .and. k1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(j,i1)
    CAP_doubly_rijk_slmn = dsqrt(1.5d0)*CAP_doubly_rijk_slmn

  endif
  
 
 elseif(typ2 == 1) then

  rt = r
  it = i
  jt = j
  kt = k

  r = r1
  i = i1
  j = j1
  k = k1

  r1 = rt
  i1 = it
  j1 = jt
  k1 = kt

  if(typ1 == 21) then

!! HERE K1 AND I ARE INTER-CHANGED IN COMPARISON OF YASEN EQUATIONS
!! I SHOULD MAKE SURE THAT IT IS OK

! old
!    if(r==r1 .and. i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - D00
!    if(i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(r,r1)
!    if(r==r1 .and. i == i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 2d0*mat_CAP_MO(j1,j)
!    if(r==r1 .and. i == i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 2d0*mat_CAP_MO(k1,k)
!    if(r==r1 .and. k == k1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(i1,i)
!    if(r==r1 .and. k == i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(k1,i)
!    if(r==r1 .and. i == j1 .and. j == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 4d0*mat_CAP_MO(i1,j)
!    if(r==r1 .and. j == i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 2d0*mat_CAP_MO(j1,i)
!    CAP_doubly_rijk_slmn = 0.5d0*CAP_doubly_rijk_slmn/dsqrt(2d0)

! new 23.02.2011
    if(r1==r .and. i1 == j .and. j1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 2d0*mat_CAP_MO(j,k1)
    if(r1==r .and. i1 == j .and. k1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0*mat_CAP_MO(j,j1)
    if(r1==r .and. j1 == i .and. k1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 2d0*mat_CAP_MO(j,i1)
    if(r1==r .and. j1 == j .and. k1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0*mat_CAP_MO(j,i1)
    if(r1==r .and. i1 == i .and. j1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0*mat_CAP_MO(j,k1)
    if(r1==r .and. i1 == i .and. k1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 1d0*mat_CAP_MO(j,j1)
    CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn/dsqrt(2d0)

  else

!! HERE K1 AND I ARE INTER-CHANGED IN COMPARISON OF YASEN EQUATIONS
!! I SHOULD MAKE SURE THAT IT IS OK

! old
!    if(r==r1 .and. i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - D00
!    if(i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(r,r1)
!    if(r==r1 .and. i == i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 2d0*mat_CAP_MO(j1,j)
!    if(r==r1 .and. i == i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 2d0*mat_CAP_MO(k1,k)
!    if(r==r1 .and. j == i1 .and. k == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(k1,i)
!    if(r==r1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(i1,i)
!    CAP_doubly_rijk_slmn = 0.5d0*dsqrt(1.5d0)*CAP_doubly_rijk_slmn

! new 23.02.2011
    if(r1==r .and. i1 == i .and. j1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(j,k1)
    if(r1==r .and. i1 == i .and. k1 == j) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(j,j1)
    if(r1==r .and. i1 == j .and. k1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(j,j1)
    if(r1==r .and. j1 == j .and. k1 == i) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(j,i1)
    CAP_doubly_rijk_slmn = dsqrt(1.5d0)*CAP_doubly_rijk_slmn

  endif

!! CHANGE BACK THE INDICES AS IT WAS IN INTENT IN

  rt = r
  it = i
  jt = j
  kt = k

  r = r1
  i = i1
  j = j1
  k = k1

  r1 = rt
  i1 = it
  j1 = jt
  k1 = kt

 elseif(typ1 == 21 .and. typ2 == 21) then

  if(r==r1 .and. i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + D00
  if(i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(r,r1)
  if(r == r1 .and. k==k1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(i1,i)
  if(r == r1 .and. k==k1 .and. i == i1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(j1,j)
  if(r == r1 .and. j==j1 .and. i == i1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(k1,k)
  if(r == r1 .and. k==k1 .and. j == i1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0*mat_CAP_MO(j1,i)
  if(r == r1 .and. k==k1 .and. i == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0*mat_CAP_MO(i1,j)
  if(r == r1 .and. j==k1 .and. i == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0*mat_CAP_MO(i1,k)
  if(r == r1 .and. k==j1 .and. j == i1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0*mat_CAP_MO(k1,i)
  if(r == r1 .and. k==j1 .and. i == i1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0*mat_CAP_MO(k1,j)
  if(r == r1 .and. j==k1 .and. i == i1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0*mat_CAP_MO(j1,k)

 elseif(typ1 == 21) then

  if(r==r1 .and. i==j1 .and. j==k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(i1,k)
  if(r==r1 .and. i==j1 .and. k==k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(i1,j)
  if(r==r1 .and. j==i1 .and. k==k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(j1,i)
  if(r==r1 .and. j==i1 .and. k==j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(k1,i)
  if(r==r1 .and. i==i1 .and. j==k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(j1,k)
  if(r==r1 .and. i==i1 .and. k==j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(k1,j)
  CAP_doubly_rijk_slmn =  0.5d0*dsqrt(3d0)*CAP_doubly_rijk_slmn

 elseif(typ2 == 21) then

  rt = r
  it = i
  jt = j
  kt = k

  r = r1
  i = i1
  j = j1
  k = k1

  r1 = rt
  i1 = it
  j1 = jt
  k1 = kt

  if(r==r1 .and. i==j1 .and. j==k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(i1,k)
  if(r==r1 .and. i==j1 .and. k==k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(i1,j)
  if(r==r1 .and. j==i1 .and. k==k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(j1,i)
  if(r==r1 .and. j==i1 .and. k==j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(k1,i)
  if(r==r1 .and. i==i1 .and. j==k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(j1,k)
  if(r==r1 .and. i==i1 .and. k==j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(k1,j)
  CAP_doubly_rijk_slmn = 0.5d0*dsqrt(3d0)*CAP_doubly_rijk_slmn

!CHANGE BACK THE INDICES AS IT WAS IN INTENT IN

  rt = r
  it = i
  jt = j
  kt = k

  r = r1
  i = i1
  j = j1
  k = k1

  r1 = rt
  i1 = it
  j1 = jt
  k1 = kt
  
 elseif(typ1 == 22 .and. typ2 == 22) then

  if(r==r1 .and. i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + D00
  if(i==i1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + mat_CAP_MO(r,r1)
  if(r==r1 .and. i == i1 .and. j == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(k1,k)
  if(r==r1 .and. i == i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(j1,j)
  if(r==r1 .and. j == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - mat_CAP_MO(i1,i)
  if(r==r1 .and. i == i1 .and. j == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 0.5d0*mat_CAP_MO(j1,k)
  if(r==r1 .and. i == j1 .and. j == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0*mat_CAP_MO(i1,k)
  if(r==r1 .and. i == i1 .and. k == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 0.5d0*mat_CAP_MO(k1,j)
  if(r==r1 .and. i == j1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 0.5d0*mat_CAP_MO(i1,j)
  if(r==r1 .and. j == i1 .and. k == k1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn - 0.5d0*mat_CAP_MO(j1,i)
  if(r==r1 .and. j == i1 .and. k == j1) CAP_doubly_rijk_slmn = CAP_doubly_rijk_slmn + 0.5d0*mat_CAP_MO(k1,i)

!  write(*,'("<",4(i1,","),"|O|",4(i1,","),">",f12.6)')i,j,k,r,i1,j1,k1,r1, CAP_doubly_rijk_slmn

 else

  write(*,*)"What is going on here, Idiot!"
  write(*,*)typ1, typ2
  stop

 endif

end function CAP_doubly_rijk_slmn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

