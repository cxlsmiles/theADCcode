program main2
use input
use CAPADC
implicit none

double complex, dimension(:,:), allocatable :: H
integer :: iCAP, nrCAP, jCAP, niCAP

double precision :: srCAP, siCAP, CAP_step, mat
double complex   :: imag
integer :: i, j

double complex, dimension(:,:), allocatable :: Hd, U, t1, t2
double complex, dimension(:), allocatable :: ed, t, t3

imag = dcmplx(0d0,1d0)

!! READS MO FROM SAJEEV CODE

read(*,*)nel
read(*,*)niCAP, CAP_step
read(*,*)srCAP


eADCmin = 0d0
eADCmax = 63.0d0/27.211d0
!eADCmax = 40.54765532000000/27.211d0
write(*,*)"READS ADC CAP MATRIX ELTS"
call read_CAP

!! Changes the strength (sCAP) of the CAP and diagonalize H = E - i*sCAP*HCAP

allocate(H(nADC,nADC))
allocate(ed(nADC),U(nADC,nADC))
allocate(t1(nADC,nADC),t2(nADC,nADC),t3(2*nADC))

write(*,*)"DIAGONALIZES CAP MATRIX"

!  CAP_step = 0.05d0
!  srCAP = 200d0

do iCAP = 1, 1  !nrCAP
  srCAP = srCAP + 0.0d0
 do jCAP = 1, niCAP !niCAP

  siCAP = (jCAP*(jCAP+1))*CAP_step
!  siCAP = jCAP*CAP_step

  H(:,:) = 0d0

  do i = 1, nADC
     H(i,i) = eADC(i)
  enddo

  do i = 1, nADC
    do j = 1, nADC
      H(i,j) = H(i,j) + srCAP*HCAP(i,j) - imag*siCAP*HCAP(i,j)
    enddo
  enddo

! call diagonalize
! DIAGONALISATION OF H USING LAPACK

   call ZGEEV('N','N', nADC, H, nADC, ed, t1, nADC, t2, nADC, U, 2*nADC, t3, i )

  do i = 1, nADC
   WRITE(25,'(3(f20.12,2X))')siCAP, ed(i)*27.211d0
  enddo

 enddo
enddo

deallocate(ed,t1,t2,t3)
deallocate(U)
deallocate(HCAP,eADC,H)

!! THE END !!

end program main2
