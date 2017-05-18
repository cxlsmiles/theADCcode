module CAPADC
use def_config
implicit none

integer :: nADC
double precision, dimension(:), allocatable :: eADC
double precision, dimension(:,:), allocatable :: HCAP
double precision :: eADCmin, eADCmax


 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_CAP
implicit none

integer :: tnADC
double precision, dimension(:), allocatable :: teADC
integer, dimension(:), allocatable :: indexit
double precision, dimension(:,:), allocatable :: tHCAP
double precision :: mat

integer :: i, j, k

open(file='capdumpfile.dat',unit=21)

read(21,*)tnADC
!write(*,*)tnADC
allocate(tHCAP(tnADC,tnADC))
allocate(teADC(tnADC), indexit(tnADC))

tHCAP(:,:) = 0d0
k = 0

do i = 1, tnADC-1
 read(21,*)teADC(i)
 if(teADC(i).ge.eADCmin .and. teADC(i) .le. eADCmax) then
   write(*,*)i,teADC(i)*27.211d0
   k = k + 1
   indexit(k) = i 
 endif
enddo

   nADC = k
   write(*,*)nADC

do i = 1, tnADC*tnADC
   read(21,*,end=100)k,j,mat
!   write(*,*)k,j,mat
   tHCAP(k+1,j+1) = mat
   tHCAP(j+1,k+1) = mat
enddo

100 close(21)

close(21)

allocate(HCAP(nADC,nADC))
allocate(eADC(tnADC))

do i = 1, nADC
  eADC(i) = teADC(indexit(i))
 do j = 1, nADC
   HCAP(i,j) = tHCAP(indexit(i),indexit(j))
   write(100,*)i,j, HCAP(i,j)
 enddo
enddo

deallocate(tHCAP,teADC,indexit)


end subroutine read_CAP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module CAPADC
