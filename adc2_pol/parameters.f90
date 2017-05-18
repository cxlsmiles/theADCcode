module parameters

  use constants

!!$*************************************
!!$*******MOLCAS provided variables*****
!!$************************************* 
!!$int*4 nIrr - number of irreps
!!$int*4 nBas - number of basis functions
!!$int*4 nCen - number of atomic centres
!!$int*4 nOcc - number of occupied orbitals
!!$int*4 nVirt - number of virtual orbitals
!!$int*4 Ehf - Hartree-Fock energy of the ground state
!!$orbSym - int*4 array(nBas) containing irrep labels of MO's
!!$e - real*8 array(nBas) containing MO's energies
!!$occNum - real*8 array(nBas) containing MO's occupation numbers 
!!$x,y,z-dipole real*8 array(nBas,nBas) containing dipole moment matrix elements 
  
  integer*4 :: nIrr,nBas,nCen
  integer :: nOcc,nVirt
  real(d) :: Ehf
  integer*4, dimension(:), allocatable :: orbSym
  real(d), dimension(:), allocatable :: e,occNum
  real(d), dimension(:,:), allocatable :: x_dipole,y_dipole,z_dipole,dpl

!!$**************************************************
!!$*********User provided variables******************
!!$**************************************************

!!$int hinit - initial hole orbital number
!!$int nirrep - number of irrep of the initial excitation
!!$int idiag - diagonalisation proc. for the init. state:1-Lapack,2-Davidson; activated if chrun is not 'direct'
!!$int fdiag - diagonalisation proc. for the final state:1-Lapack,2-Lanczos; activated if chrun is not 'direct' 
!!$int method - comp. method: (0) - fano, (1) - tda, (2) - adc2 (3) - adc2e
!!$int fmethod - fano method: (1) - tda, (2) - adc2 (3) - adc2e
!!$int davstates - number of requested Davidson initial states 
!!$int lancstates - number of lanczos computed final states; Number of printed states in simple ADC2 calculations
!!$int numinista - number of initial states in Fano calculations;
!!$int array(1:numinista) ninista - the initial states for which appropriate Lanczos generated final state manifold 
!!$is to be computed  
!!$int stiprilev - print level in the Stieltjes subroutine
!!$int norder - Stieltjes order
!!$real*8 minc - minimum value of accepted adc matr. elem.
!!$real*8 mspacewi - minimum weight of single exc. in an allowed initial state
!!$real*8 mspacewf - minimum weight of single exc. in an allowed final state
!!$logical readband - a flag requesting the final eigenvectors in an energy band (eupper,elower)
!!$real*8 eupper,elower - energy band borders 
!!$roccnum - int array(nBas) rearranged MO's in the order occ -> virt.
!!$int array(8,8) MT -  multiplication table for the abelian groups d2h, etc.
!!$int array(nhcentre) hcentre - contains labels of the occupied orbitals on an atom carrying the initial hole
!!$char(4) chrun - 'direct'-Lapack diagonalization of both matrices, 'save'-save matrices 
!!$needed for either Lanczos or Davidson routine, 'read' read vectors produced by the external 
!!$diag. routine and get lifetimes 
!!$char(36) davname,lancname  - the names of the davidson and lanczos vector files
!!$char(1) tranmom  - dipole moment index 'x','y', or 'z' corresp. to the nirrep.
!!$integer array(1:lmain) stvc_lbl  - Damit wirder der Lanc-Startblock festgelegt
!!$integer info if 1 stops execution after printing the configuration tables
!!$integer ninista gives the number of the fanostate among davidson eigenvectors
  character(1) :: tranmom,tranflag
  character(4) :: chrun
  character(36) :: lancname,davname
  integer :: hinit,nirrep,idiag,fdiag,method,fmethod,davstates,lancstates,stiprilev,numinista,norder,info
  integer, parameter :: nhcentre=20
  integer, dimension(0:nhcentre) ::  hcentre
  integer, parameter :: nhneighb=20
  integer, dimension(0:nhneighb) ::  hneighb
  real(d) :: minc,mspacewi,mspacewf,eupper,elower
  integer, dimension(:), pointer :: roccnum
  integer :: mgvdim
  real(d), dimension(:), allocatable :: mgvec
  integer, dimension(8,8) :: MT
  logical :: readband
  integer, dimension(1000) :: stvc_lbl
  integer :: ninista

!!$************************************************
!!$**********Physical Cobnstants*******************
!!$************************************************

!!$ abohr - Bohr radius [cm]
!!$ fsconstinv - inverted fine structure constant
!!$ os2cs - oscillator strength [a.u.] to cross-section [Mb] conversion factor
!!$ omega - photon energy, required by Stieltjes_phi, photoionisation routine.

  real(d), parameter :: abohr=5.2918e-9
  real(d), parameter :: fsconstinv=137._d
  real(d), parameter :: os2cs=4.0347443
  real(d), parameter :: omega=3.0_d
  
!!$************************************************
!!$**********Lanczos Parameters********************
!!$************************************************  

  character*4 :: mtxidl
  integer :: ncycles,maxmem,memx,lmain
  integer :: mode,nprint
  integer :: maxiter
  integer, dimension(2) :: iparm
  real(d) :: wthr
  real(d), dimension(2) :: erange
  real(d) :: unit
  real(d),dimension(5) :: fparm

!!$************************************************
!!$**********Davidson Parameters********************
!!$************************************************  

  character*4 :: mtxidd
  integer :: dmain
  logical :: myb0,transp




end module parameters
