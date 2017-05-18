module ADC
use def_config
use input
implicit none

type ADC_eigenvectors

 integer :: nconfig
 double precision :: eigval
 double precision, dimension(:), allocatable :: coeff
 type(spin_adapted_config), dimension(:), allocatable :: cfg

end type ADC_eigenvectors

integer :: nADC_eigvec
integer, dimension(:), allocatable :: indexADC
integer :: nADCselect
double precision, dimension(:,:), allocatable :: HCAPADC

 contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function convert(a,length)
implicit none

character(*), intent(in) :: a
integer, intent(in) :: length
integer :: convert

integer :: i

 convert = 0
 do i = 1, length
  convert = convert + 10**(length-i)*(ichar(a(i:i))-48)
 enddo 

end function convert


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_eigvec(eigvec)
implicit none

type(ADC_eigenvectors), dimension(:), intent(inout) :: eigvec

type charac
  character(20) :: c
end type charac

type(charac), dimension(10) :: a, b

integer :: i, j, k
integer :: length!, convert

nADCselect = 0

open(file='capdumpfile.dat',unit=10)

  read(10,*) i
  write(*,*)i

do i = 1, nADC_eigvec

  read(10,*)eigvec(i)%eigval
  if(eigvec(i)%eigval .ge. eADCmin .and. eigvec(i)%eigval .le. eADCmax) then
    nADCselect = nADCselect + 1
    indexADC(nADCselect) = i
  endif
  read(10,*)eigvec(i)%nconfig
!    write(*,*)i,eigvec(i)%eigval,eigvec(i)%nconfig
  if(.not.allocated(eigvec(i)%coeff)) allocate(eigvec(i)%coeff(eigvec(i)%nconfig))
  if(.not.allocated(eigvec(i)%cfg)) allocate(eigvec(i)%cfg(eigvec(i)%nconfig))

  do j = 1, eigvec(i)%nconfig
    read(10,*)eigvec(i)%cfg(j)%nholes,eigvec(i)%cfg(j)%nparts,(a(k)%c,k=1,eigvec(i)%cfg(j)%nholes), &
    (b(k)%c,k=1,eigvec(i)%cfg(j)%nparts),eigvec(i)%cfg(j)%typ,eigvec(i)%coeff(j)


! must be changed back later with proper Yasen output
!    read(*,*)eigvec(i)%cfg(j)%nholes,eigvec(i)%cfg(j)%nparts,(a(k)%c,k=1,eigvec(i)%cfg(j)%nholes-1), &
!    (b(k)%c,k=1,eigvec(i)%cfg(j)%nparts),eigvec(i)%cfg(j)%typ,eigvec(i)%coeff(j)
!
!    eigvec(i)%cfg(j)%nholes = eigvec(i)%cfg(j)%nholes - 1
!!!


     if(eigvec(i)%cfg(j)%nholes==1)  then

       length=index(a(1)%c,' ')-1
       eigvec(i)%cfg(j)%ihole(1) = ichar(a(1)%c(2:length))-48

!       write(*,'(3(i2,1x),4x,a1,1(i3,1x),a1,2x,i2,4x,f12.6)')j-1,eigvec(i)%cfg(j)%nholes,eigvec(i)%cfg(j)%nparts,"<",(eigvec(i)%cfg(j)%ihole(k),k=1,eigvec(i)%cfg(j)%nholes), &
!       (eigvec(i)%cfg(j)%ipart(k),k=1,eigvec(i)%cfg(j)%nparts),"|", eigvec(i)%cfg(j)%typ, eigvec(i)%coeff(j)

     endif


     if(eigvec(i)%cfg(j)%nholes==2)  then

       length=index(a(1)%c,' ')-1
       eigvec(i)%cfg(j)%ihole(1) = ichar(a(1)%c(2:length))-48

       length=index(a(2)%c,' ')-1
       eigvec(i)%cfg(j)%ihole(2) = convert(a(2)%c,length)

       length=index(b(1)%c,'|')-1
       eigvec(i)%cfg(j)%ipart(1) = convert(b(1)%c(:),length)

!       write(*,'(3(i2,1x),4x,a1,3(i3,1x),a1,2x,i2,4x,f12.6)')j-1, eigvec(i)%cfg(j)%nholes,eigvec(i)%cfg(j)%nparts,"<",(eigvec(i)%cfg(j)%ihole(k),k=1,eigvec(i)%cfg(j)%nholes), &
!       (eigvec(i)%cfg(j)%ipart(k),k=1,eigvec(i)%cfg(j)%nparts),"|", eigvec(i)%cfg(j)%typ, eigvec(i)%coeff(j)

     endif

  enddo

! stop
enddo


close(10)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_eigvec_doubly(eigvec)
implicit none

type(ADC_eigenvectors), dimension(:), intent(inout) :: eigvec

type charac
  character(20) :: c
end type charac

type(charac), dimension(10) :: a, b

integer :: i, j, k
integer :: length!, convert
double precision, pointer :: ptr

open(file='capdumpfile.dat',unit=10)

  read(10,*) i
  write(*,*)i

do i = 1, nADC_eigvec

  read(10,*)eigvec(i)%eigval
  write(*,*)eigvec(i)%eigval
  if(eigvec(i)%eigval .ge. eADCmin .and. eigvec(i)%eigval .le. eADCmax) then
    nADCselect = nADCselect + 1
    indexADC(nADCselect) = i
  endif
  read(10,*)eigvec(i)%nconfig

  if(.not.allocated(eigvec(i)%coeff)) allocate(eigvec(i)%coeff(eigvec(i)%nconfig))
  if(.not.allocated(eigvec(i)%cfg)) allocate(eigvec(i)%cfg(eigvec(i)%nconfig))


  do j = 1, eigvec(i)%nconfig
    read(10,*)eigvec(i)%cfg(j)%nholes,eigvec(i)%cfg(j)%nparts,(a(k)%c,k=1,eigvec(i)%cfg(j)%nholes), &
    (b(k)%c,k=1,eigvec(i)%cfg(j)%nparts),eigvec(i)%cfg(j)%typ,eigvec(i)%coeff(j)

     if(eigvec(i)%cfg(j)%nholes==2)  then

       length=index(a(1)%c,' ')-1
       eigvec(i)%cfg(j)%ihole(1) = ichar(a(1)%c(2:length))-48

       length=index(a(2)%c,'|')-1
       eigvec(i)%cfg(j)%ihole(2) = convert(a(2)%c,length)
!       eigvec(i)%cfg(j)%ihole(2) = ichar(a(2)%c(1:length))-48

!       write(*,'(3(i2,1x),4x,a1,2(i3,1x),a1,2x,i2,4x,f12.6)')j-1,eigvec(i)%cfg(j)%nholes,eigvec(i)%cfg(j)%nparts,"<",(eigvec(i)%cfg(j)%ihole(k),k=1,eigvec(i)%cfg(j)%nholes), &
!       (eigvec(i)%cfg(j)%ipart(k),k=1,eigvec(i)%cfg(j)%nparts),"|", eigvec(i)%cfg(j)%typ, eigvec(i)%coeff(j)

     endif


     if(eigvec(i)%cfg(j)%nholes==3)  then

       length=index(a(1)%c,' ')-1
       eigvec(i)%cfg(j)%ihole(1) = ichar(a(1)%c(2:length))-48

       length=index(a(2)%c,' ')-1
       eigvec(i)%cfg(j)%ihole(2) = convert(a(2)%c,length)
!       eigvec(i)%cfg(j)%ihole(2) = ichar(a(2)%c(1:length))-48

       length=index(a(3)%c,' ')-1
       eigvec(i)%cfg(j)%ihole(3) = convert(a(3)%c,length)
!       eigvec(i)%cfg(j)%ihole(3) = ichar(a(3)%c(1:length))-48

       length=index(b(1)%c,'|')-1
       eigvec(i)%cfg(j)%ipart(1) = convert(b(1)%c(:),length)

!       write(*,'(3(i2,1x),4x,a1,4(i3,1x),a1,2x,i2,4x,f12.6)')j-1, eigvec(i)%cfg(j)%nholes,eigvec(i)%cfg(j)%nparts,"<",(eigvec(i)%cfg(j)%ihole(k),k=1,eigvec(i)%cfg(j)%nholes), &
!       (eigvec(i)%cfg(j)%ipart(k),k=1,eigvec(i)%cfg(j)%nparts),"|", eigvec(i)%cfg(j)%typ, eigvec(i)%coeff(j)

     endif

  enddo

! stop
enddo

 close(10)
WRITE(*,*)"END ADC EIGENVECTORS"

end subroutine read_eigvec_doubly

end module ADC
