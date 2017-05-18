!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Computes the CAP between two ADC eigenvectors

subroutine ADCCAP_SINGLY(eigvec1,eigvec2)
use ADC
use def_config
use CAP_MO
implicit none

 type(ADC_eigenvectors), intent(in) :: eigvec1, eigvec2

 double precision :: a

 double precision :: mat_spac
 integer :: i, j

 do i = 1, eigvec1%nconfig
   do j = 1, eigvec2%nconfig

       call CAP_spin_adapted_config(eigvec1%cfg(i),eigvec2%cfg(j),mat_spac)
       HCAPADC(i,j) = mat_spac

   enddo
 enddo

end subroutine ADCCAP_SINGLY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Computes the CAP between two ADC eigenvectors

subroutine CAP_SINGLY(eigvec1,eigvec2,mat)
use ADC
use def_config
use CAP_MO
implicit none

 type(ADC_eigenvectors), intent(in) :: eigvec1, eigvec2
 double precision, intent(out) :: mat

 double precision :: a

 double precision :: mat_spac
 integer :: i, j

 mat = 0d0

 do i = 1, eigvec1%nconfig
   do j = 1, eigvec2%nconfig

      mat = mat + eigvec1%coeff(i)*eigvec2%coeff(j)*HCAPADC(i,j)

   enddo
 enddo

end subroutine CAP_SINGLY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Computes the CAP between two spin adapted config.

subroutine CAP_spin_adapted_config(spac1,spac2,mat_spac)
use def_config
implicit none

 type(spin_adapted_config), intent(inout) :: spac1, spac2
 double precision, intent(out) :: mat_spac

 double precision :: CAP_i_j, CAP_i_rkl, CAP_rij_slm
 integer :: i, j

! only doublet spin

 if(spac1%nholes==1 .and. spac2%nholes==1) then

    mat_spac = CAP_i_j(spac1%ihole(1),spac2%ihole(1))

 elseif(spac1%nholes==1) then

    mat_spac = CAP_i_rkl(spac1%ihole(1),spac2%ipart(1),spac2%ihole(1),spac2%ihole(2),spac2%typ)

 elseif(spac2%nholes==1) then

    mat_spac = CAP_i_rkl(spac2%ihole(1),spac1%ipart(1),spac1%ihole(1),spac1%ihole(2),spac1%typ)

 else

    mat_spac = CAP_rij_slm(spac1%ipart(1),spac1%ihole(1),spac1%ihole(2),&
      spac1%typ,spac2%ipart(1),spac2%ihole(1),spac2%ihole(2),spac2%typ)

 endif

end subroutine CAP_spin_adapted_config

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_i_j(i,j)
use CAP_MO
implicit none

 double precision :: CAP_i_j
 integer, intent(in) :: i, j

  CAP_i_j = 0d0

  if(i==j) CAP_i_j = D00

  CAP_i_j =  CAP_i_j - mat_CAP_MO(i,j)

end function CAP_i_j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_i_rkl(i,r,k,l,typ)
use CAP_MO
implicit none

 double precision :: CAP_i_rkl
 integer, intent(in) :: i, r, k, l, typ

 CAP_i_rkl = 0d0


 if(typ==1) then

 ! k = l

   if(i==k) CAP_i_rkl = CAP_i_rkl + mat_CAP_MO(k,r)

 elseif(typ==21) then

   if(i==l) CAP_i_rkl = CAP_i_rkl + mat_CAP_MO(k,r)/dsqrt(2d0)
   if(i==k) CAP_i_rkl = CAP_i_rkl + mat_CAP_MO(l,r)/dsqrt(2d0)

 elseif(typ==22) then

   if(i==l) CAP_i_rkl = CAP_i_rkl + mat_CAP_MO(k,r)*dsqrt(1.5d0)
   if(i==k) CAP_i_rkl = CAP_i_rkl - mat_CAP_MO(l,r)*dsqrt(1.5d0)

 endif

end function CAP_i_rkl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CAP_rij_slm(b,i,j,typ1,a,k,l,typ2)
use CAP_MO
implicit none

 double precision :: CAP_rij_slm
 integer :: a, b, i, j, k, l, typ1, typ2
 integer :: a1, b1, i1, j1, k1, l1
 integer :: at, bt, it, jt, kt, lt

 CAP_rij_slm = 0d0

 if(typ1==1 .and. typ2==1) then

  if(b==a .and. i==k) CAP_rij_slm = CAP_rij_slm + D00 - 2d0*mat_CAP_MO(i,k)
  if(i==k) CAP_rij_slm = CAP_rij_slm + mat_CAP_MO(a,b)

 elseif(typ1==1 .and. typ2==21) then 

! i = j

!  if(a==b .and. i==k .and. i==l) CAP_rij_slm = CAP_rij_slm + D00/dsqrt(2d0)
!  if(i==k .and. i==l) CAP_rij_slm = CAP_rij_slm + mat_CAP_MO(b,a)/dsqrt(2d0)
  if(b==a .and. i==k) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(i,l)
  if(b==a .and. i==l) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(i,k)

  CAP_rij_slm = CAP_rij_slm*dsqrt(2d0)

 elseif(typ1==21 .and. typ2==1) then

  at = a
  it = i
  jt = j

  a = b
  i = k
  j = l

  b = at
  k = it
  l = jt

! k = l

!  if(a==b .and. i==k .and. i==l) CAP_rij_slm = CAP_rij_slm + D00/dsqrt(2d0)
!  if(i==k .and. i==l) CAP_rij_slm = CAP_rij_slm + mat_CAP_MO(b,a)/dsqrt(2d0)
  if(b==a .and. i==k) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(i,l)
  if(b==a .and. i==l) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(i,k)

  CAP_rij_slm = CAP_rij_slm*dsqrt(2d0)
  
  at = a
  it = i
  jt = j

  a = b
  i = k
  j = l

  b = at
  k = it
  l = jt

 elseif(typ1==1 .and. typ2==22) then 

  CAP_rij_slm = 0d0

 elseif(typ1==22 .and. typ2==1) then

  CAP_rij_slm = 0d0

 elseif(typ1==21 .and. typ2==22) then 

  CAP_rij_slm = 0d0

 elseif(typ1==22 .and. typ2==21) then

  CAP_rij_slm = 0d0

 elseif(typ1==21 .and. typ2==21) then

  if(a==b .and. i==k .and. j==l) CAP_rij_slm = CAP_rij_slm + D00 
  if(i==k .and. j==l) CAP_rij_slm = CAP_rij_slm + mat_CAP_MO(b,a)
  if(a==b .and. i==k) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(j,l)
  if(a==b .and. i==l) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(k,j)
  if(a==b .and. j==l) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(k,i)
  if(a==b .and. k==j) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(l,i)


 elseif(typ1==22 .and. typ2==22) then

  if(a==b .and. i==k .and. j==l) CAP_rij_slm = CAP_rij_slm + D00 
  if(i==k .and. j==l) CAP_rij_slm = CAP_rij_slm + mat_CAP_MO(b,a)
  if(a==b .and. i==k) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(j,l)
  if(a==b .and. i==l) CAP_rij_slm = CAP_rij_slm + mat_CAP_MO(k,j)
  if(a==b .and. j==l) CAP_rij_slm = CAP_rij_slm - mat_CAP_MO(k,i)
  if(a==b .and. k==j) CAP_rij_slm = CAP_rij_slm + mat_CAP_MO(l,i)

 endif

end function CAP_rij_slm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

