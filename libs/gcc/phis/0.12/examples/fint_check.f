      Program FInt_Check
      
      Implicit None

      Integer,Parameter :: backend=131072
      Real(8),Parameter :: zero=0.D0,half=0.5D0
      Integer           :: cap,pos,p,q,r,s=1,nsym,nbas,natoms,i,j
      Real(8)           :: v_pqrs,iiii,jjjj,iijj,ijij

!...External Functions...
      Integer, External :: phis_init
      Real(8), External :: vpqrs

!Indirect way for extracting 2-electron integrals

      write(6,'(A)') "   p   q   r   s |    Vpqrs"
      write(6,'(A)') "-----------------+-------------"
!Loading all 2-electron integrals into RAM
      cap=phis_init(backend)
      call phis_load_vpqrs    
      call phis_get_info(nsym,nbas,natoms)
!      call phis_mc_epsi .... !check SPEC file
!Checking integrals according to Roothan`s relation (26):
! C.C.J. Roothan "New Developments in Molecular Orbital Theory"  //
! Reviews of modern physics, V. 23, N. 2, P. 69-89, P. 72
      do i=1,nbas
         do j=1,i,nbas
            iiii=vpqrs(i,i,i,i) ! J[ii]
            jjjj=vpqrs(j,j,j,j) ! J[jj]
            iijj=vpqrs(i,i,j,j) ! J[ij]
            ijij=vpqrs(i,j,i,j) ! K[ij]
            if (.not.((ijij<=iijj).and.(iijj<=half*(iiii+jjjj)))) Stop "Something is wrong"
         end do
      end do 
         
      End Program FInt_Check
