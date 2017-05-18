module main_cap
use input
use ADC
use def_config
use CAP_MO
use test
implicit none

character(80) :: fileout

type(ADC_eigenvectors), dimension(:), allocatable :: eigvec
double precision, dimension(:,:), allocatable :: HCAP
double complex, dimension(:,:), allocatable :: H
integer :: iCAP, nrCAP, jCAP, niCAP, offsetCAP

double precision :: srCAP, siCAP, CAP_step, mat
double complex   :: imag
integer :: i, j, k

double complex, dimension(:,:), allocatable :: Hd, U, t1, t2
double complex, dimension(:), allocatable :: ed, t, t3

character(3) :: incr

! FROM AND FOR SAJEEV CODE

INTEGER *4 ISPC, NB
PARAMETER (ISPC=3000)

contains
  
  subroutine docap(IO,dens,engs,dipMO)
    
    TYPE INFO
       INTEGER :: NMO,NOCC, NICAP, IOFFSET, INCR, NHOLES, MULT
       DOUBLE PRECISION :: SRCAP, SICAP, EADCMAX, EADCMIN
       CHARACTER*64 OUTPUT
    END TYPE INFO
    
    TYPE (INFO), intent(in) :: IO
  
    integer :: nmos, nel
    double precision, intent(in) :: dens(IO%nmo,IO%nmo)
    double precision, intent(in) :: dipMO(IO%nmo,IO%nmo)
    double precision, intent(in) :: engs(IO%nmo)
    double precision :: amat(IO%nmo,IO%nmo)


    nmos = IO%nmo
    nmo  = IO%nmo
    allocate(DMO(nmos,nmos),EMO(nmos))
    EMO(:) = engs(:)
    DMO(:,:) = dens(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! READ INPUT PARAMETERS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! The input paramaters are passed by the C++ code
    !    read(*,*)nel ! number of active electrons
    !    read(*,*)niCAP, offsetCAP, CAP_step, incr ! number of trajectories, iteration offset, step, exponential or linear increase
    !    read(*,*)srCAP ! strength of the real part of the CAP
    niCAP = IO%nicap
    offsetCAP = IO%ioffset
    CAP_step = IO%sicap
    srCAP = IO%srcap
    select case (IO%incr)
    case (1)
       incr = 'exp'
    case (2)
       incr = 'lin'
    end select

    nocc = IO%nocc
    nel = 2*nocc
    nvirt = nMO - nocc
    
    !read(*,*)eADCmin, eADCmax ! Limit the energy range (in eV) of the ADC eigenvectors used
    !read(*,*)fileout ! name of the file where trajectories will be printed
    eADCmin = IO%eadcmin
    eADCmax = IO%eadcmax
    fileout = IO%output


    !! READS MO FROM SAJEEV CODE 
    ! Done now in the c++   part
    allocate(mat_CAP_MO(nMOs,nMOs))
    mat_CAP_MO(:,:) = dipMO(:,:)

! WE DO NOT WANT TO APPLY CAP TO THE OCC MO

 do i = 1, nel/2
    mat_CAP_MO(i,:) = 0d0
    mat_CAP_MO(:,i) = 0d0
 enddo

    ! OLD:
    ! All we need here is dipMO
    ! READS SOME INFOS FROM GUS OUTPUT FILES
    !    open(unit=5,file='capinput.txt',STATUS='UNKNOWN',FORM='FORMATTED')
    !    open(unit=8,file='mocoef.txt',STATUS='UNKNOWN',FORM='FORMATTED')
    !    read (8, *) NB
    !    close(8)
    !COMPUTES CAP MATRIX ELEMENTS IN MO BASIS
    !AND PRINTS THEM IN capmat.dat
    !    call Compute_CAP_MO(NB,ISPC)
    !    call read_CAP_MO
    !     allocate(mat_CAP_MO(0:nMOs,0:nMOs))
    !     nMO = nMOs
    !     mat_CAP_MO(:,:) = transpose(dipMO(0:nMOs,0:nMOs))
         open(unit=10,file='MOdip_mom.dat')
         do i = 1, nMOs
             write(10,'(6(2(i4,1X),1f12.8,1X))')(i,j,mat_CAP_MO(i,j),j=1,nMOs)
         enddo
         close(10)
    !    write (*,*) dipMO(:,:)
    !    stop


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !NOW WE CAN START TRAJECTORIES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    imag = dcmplx(0d0,1d0)
    
      D00=0d0
!at zeroth order
    do i = 1, nel/2
       D00 = D00 + 2d0*mat_CAP_MO(i,i)
    enddo
!at second order
      D02=0d0
    do i = 1, nMO
      do j = 1, nMO
         D02 = D02 + 2d0*dens(i,j)*mat_CAP_MO(j,i)
      enddo
    enddo
! test
!    D02 = D00

    !! READS ADC EIGENVECTORS FROM YASEN CODE
    
    open(file='capdumpfile.dat',unit=10)
    read(10,*) nADC_eigvec
    write(*,*) nADC_eigvec
    close(10)
    
    allocate(eigvec(nADC_eigvec))
    allocate(indexADC(nADC_eigvec))
    
    write(*,*)"READS ADC EIGENVECTORS"
    if (IO%NHOLES .eq. 1) then 
       write(*,*)"SINGLE IONIZATION"
       call read_eigvec(eigvec)
    else 
       write(*,*)"DOUBLE IONIZATION"
       call read_eigvec_doubly(eigvec) 
    endif
    
    !! Print selected ADC eigenvectors and their main contribution
    
!    do i = 1, nADCselect
!       j = maxloc(abs(eigvec(indexADC(i))%coeff(:)),1)
!       if(eigvec(i)%cfg(j)%nholes==1)  then
!          write(*,'(i4,3X,2(f12.6,1x),4x,a1,1(i3,1x),a1)') i, eigvec(indexADC(i))%eigval,maxval(abs(eigvec(indexADC(i))%coeff(:))),"<",(eigvec(indexADC(i))%cfg(j)%ihole(k),k=1,eigvec(indexADC(i))%cfg(j)%nholes), &
!               (eigvec(indexADC(i))%cfg(j)%ipart(k),k=1,eigvec(indexADC(i))%cfg(j)%nparts),"|"
!       else
!          write(*,'(i4,3X,2(f12.6,1x),4x,a1,3(i3,1x),a1)') i, eigvec(indexADC(i))%eigval,maxval(abs(eigvec(indexADC(i))%coeff(:))),"<",(eigvec(indexADC(i))%cfg(j)%ihole(k),k=1,eigvec(indexADC(i))%cfg(j)%nholes), &
!               (eigvec(indexADC(i))%cfg(j)%ipart(k),k=1,eigvec(indexADC(i))%cfg(j)%nparts),"|"
!       endif
!    enddo
    
    allocate(HCAP(nADC_eigvec,nADC_eigvec))
    
    !! Computes CAP in ADC eigenvectors basis
    
    HCAP(:,:) = 0d0
    
    write(*,*)"COMPUTES CAP MATRIX ELEMENTS"
    if (IO%NHOLES .eq. 2) then 
    if (IO%MULT .eq. 3)  then
         call test_isr_print_triplet(eigvec(indexADC(1)),eigvec(indexADC(1)),mat)
    elseif (IO%MULT .eq. 1)  then
         call test_isr_print(eigvec(indexADC(1)),eigvec(indexADC(1)),mat)
    else
         write(*,*)"don't know this spin multiplicity,I stop",IO%MULT
         stop
    endif
    endif
          
    allocate(HCAPADC(eigvec(indexADC(1))%nconfig,eigvec(indexADC(1))%nconfig))

    if (IO%NHOLES .eq. 1) then
       write(*,*)"SINGLE IONIZATION NOT IMPLEMENTED YET, STOP"
          call ADCCAP_SINGLY(eigvec(indexADC(1)),eigvec(indexADC(1)))
!       stop
    else
       write(*,*)"DOUBLE IONIZATION"
      if (IO%MULT .eq. 1)  then
          call ADCCAP_DOUBLY(eigvec(indexADC(1)),eigvec(indexADC(1)))
      elseif (IO%MULT .eq. 3) then
          call ADCCAP_DOUBLY_TRIPLET(eigvec(indexADC(1)),eigvec(indexADC(1)))
      endif
    endif


!    write(*,*)nADCselect
!    do i = 1, nADC_eigvec
!      do j = 1, nADC_eigvec
!    open(unit=25,file=fileout)

  do i = 1, nADCselect
       do j = 1, nADCselect

    if (IO%NHOLES .eq. 1) then
!       write(*,*)"SINGLE IONIZATION"
       call CAP_SINGLY(eigvec(indexADC(i)),eigvec(indexADC(j)),mat)
    else
!       write(*,*)"DOUBLE IONIZATION",IO%MULT 
       if (IO%MULT .eq. 1) then
           call CAP_doubly(eigvec(indexADC(i)),eigvec(indexADC(j)),mat)
       elseif (IO%MULT .eq. 3) then
           call CAP_doubly_triplet(eigvec(indexADC(i)),eigvec(indexADC(j)),mat)
       endif
    endif

       HCAP(i,j) = mat
       write(25,*)i,j, HCAP(i,j),(-HCAP(i,j)+0.975)*2.54
!      stop

       enddo
    enddo

!     close(25)
!
!    STOP
    
    !! Changes the strength (sCAP) of the CAP and diagonalize H = E - i*sCAP*HCAP
    
    !allocate(H(nADC_eigvec,nADC_eigvec))
    !allocate(ed(nADC_eigvec),U(nADC_eigvec,nADC_eigvec))
    !allocate(t1(nADC_eigvec,nADC_eigvec),t2(nADC_eigvec,nADC_eigvec),t3(2*nADC_eigvec))
    
    allocate(H(nADCselect,nADCselect))
    allocate(ed(nADCselect),U(nADCselect,nADCselect))
    allocate(t1(nADCselect,nADCselect),t2(nADCselect,nADCselect),t3(2*nADCselect))
    
    write(*,*)"DIAGONALIZES CAP MATRIX"
    
    open(unit=25,file=fileout)
    
    do iCAP = 1, 1  !nrCAP
       srCAP = srCAP + 0.0d0
       do jCAP = offsetCAP, niCAP+offsetCAP !niCAP
          
          if(incr=='exp') siCAP = (jCAP*(jCAP+1))*CAP_step
          if(incr=='lin') siCAP = jCAP*CAP_step
          
          H(:,:) = 0d0
          
          !  do i = 1, nADC_eigvec
          do i = 1, nADCselect
             H(i,i) = eigvec(indexADC(i))%eigval/27.211d0
          enddo
          
          !  do i = 1, nADC_eigvec
          !    do j = 1, nADC_eigvec
          do i = 1, nADCselect
             do j = 1, nADCselect
                H(i,j) = H(i,j) + srCAP*HCAP(i,j) - imag*siCAP*HCAP(i,j)
             enddo
          enddo
          
          ! call diagonalize
          ! DIAGONALISATION OF H USING LAPACK
          
          !   call ZGEEV('N','N', nADC_eigvec, H, nADC_eigvec, ed, t1, nADC_eigvec, t2, nADC_eigvec, U, 2*nADC_eigvec, t3, i )
          
          !   call ZGEEV('N','N', nADCselect, H, nADCselect, ed, t1, nADCselect, t2, nADCselect, U, 2*nADCselect, t3, i )
          
          call ZGEEV('N','V', nADCselect, H, nADCselect, ed, t1, nADCselect, t2, nADCselect, U, 2*nADCselect, t3, i )
          U(:,:) = dconjg(t2(:,:))
          
          do i = 1, nADCselect
             mat = 0d0
             do j = 1, nADCselect
                do k = 1, nADCselect
                   mat = mat + U(j,i)*U(k,i)*srCAP*HCAP(j,k)
                enddo
             enddo
             
!             ed(i) = ed(i) - mat
          enddo
          
          
          !  do i = 1, nADC_eigvec
          do i = 1, nADCselect
             !WRITE(25,'(3(f40.12,2X))')siCAP, ed(i)*27.211d0
             WRITE(25,'(2(f40.12,2X))') ed(i)*27.211d0
          enddo
          
       enddo
    enddo
    
    close(25)
    
    deallocate(ed,t1,t2,t3)
    deallocate(U)
    deallocate(eigvec)
    deallocate(indexADC)
    deallocate(HCAP,H)
    deallocate(EMO,DMO)
    
    !! THE END !!
    
  end subroutine docap
  
end module main_cap
