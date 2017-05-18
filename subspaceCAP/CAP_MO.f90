module CAP_MO
use input
implicit none




integer :: nMO
double precision, dimension(:,:), allocatable :: mat_CAP_MO
double precision, dimension(:), allocatable :: EMO
double precision, dimension(:,:), allocatable :: DMO
! D00 =  HF CAP matrix element
double precision :: D00
double precision :: D02

 contains

subroutine read_MO

 integer :: i, j
 character(40) :: a

 open(unit=10,file='orbital.txt')
 
 do i = 1, 7
    read(10,*)a
 enddo

 do i = 1, nMO
  read(10,*)j,a,EMO(i)
  write(*,*)j,a,EMO(i)
 enddo

 close(10)

end subroutine read_MO


subroutine read_CAP_MO

 integer :: i, j

! open(unit=10,file='capmat.dat')
 open(unit=10,file='MOdip_mom.dat')

  read(10,*)nMO

  allocate(mat_CAP_MO(0:nMO,0:nMO))

  mat_CAP_MO(:,:) = 0d0

!! MAKE SURE THE ORDER OF MO ARE THE SAME IN CAP MO SAJEEV CALC.
!! AND IN YASEN CODE.

  do i = 1, nMO
    do j = 1, i
       read(10,*)mat_CAP_MO(i,j)
       mat_CAP_MO(j,i) = mat_CAP_MO(i,j)
    enddo
  enddo

 close(10)

  do i = 1, nMO
       write(10,'(6(2(i4,1X),1f20.16,1X))')(i,j,mat_CAP_MO(i,j),j=1,nMO)
  enddo

! WE DO NOT WANT TO APPLY CAP TO THE OCC MO

! do i = 1, nel/2
!    mat_CAP_MO(i,:) = 0d0
!!    mat_CAP_MO(:,i) = 0d0
! enddo

! HF CAP matrix element

 do i = 1, nel/2
    D00 = D00 + 2d0*mat_CAP_MO(i,i)
 enddo


end subroutine read_CAP_MO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deallocate_CAP_MO
deallocate(mat_CAP_MO)
end subroutine deallocate_CAP_MO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module CAP_MO
