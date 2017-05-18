!     fdump_all.f
!
!     Small sample application, which dumps all the available
!     information of the SCF calculation to stdout (6).
!     
!     JB(3/02)

      program fdump_all

      implicit none

      integer maxsym,maxbas,maxato,max_cc
      integer tmpbas,tmp_cc
      character(len=80) path
      character(len=0)  NULL
      parameter(maxsym=8,maxbas=1024,maxato=100,max_cc=13)

      Integer nBasSum, nActSum, nAtom
      Integer MapToExternal(5000)
      Integer myBas(8), nOrb(8), nOcc(8)
      Integer nAct(8), nFro(8), nDel(8)

      Real*8 Coor(3,500),Znuc(500),Escf

      integer cap,nsym,nbas,natoms

      integer poly(3,maxbas)
      integer nmb_cc(maxbas)
      real*8 cc(max_cc,maxbas)
      real*8 alpha(max_cc,maxbas)
      integer center(maxbas)

      integer i,j

!...EXTERNAL FUNCTIONS...

      integer  phis_init
      external phis_init

      Common /PHIS_MC/ nBasSum, nActSum,
     &     Escf,
     &     myBas, nOrb, nOcc, nAct, nFro, nDel,
     &     MapToExternal,
     &     Coor, Znuc, nAtom


!...FORMATS...

 100  format (i4,' |',i4,'  |',3i2,' | ',i4,'   |',3f8.3,' |',3f7.3)

! FIXME: just guk implemented
      path="/tmpa/vitja/GUK"
!!      cap=phis_init(1,trim(path))

      cap=phis_init(131072)
      call phis_get_info(nsym,nbas,natoms)
      print *, "escf",escf
      print *, "nBas",myBas
      print *, "nOrb",nOrb
      print *,"cap=",cap
      
      write(6,'(a,i4)')  'no. of irreps:         ',nsym
      write(6,'(a,i4)')  'no. of basis functions:',nbas
      write(6,'(a,i4/)') 'no. of atoms:          ',natoms
      
!----------------------------------------------------------------------C

      tmpbas=maxbas
      tmp_cc=max_cc
      call phis_get_ao(poly,nmb_cc,cc,alpha,center,tmpbas,tmp_cc)

      if(tmpbas.le.0 .or. tmp_cc.le.0) then
         write(6,'(a/a,i5/a,i5)') 'arrays too small:',     
     &         '  max_cc must be greater than',abs(tmp_cc),
     &         '  maxato must be greater than',abs(tmpbas)
         stop
      endif

      write(6,'(a/a)') '  AO |center| i j k | nmb_cc | first three '// 
     &      'alphas      | first three ccs','-----+------+-------+--'// 
     &      '------+-------------------------+----------------------'
      
      do i=1,nbas

!     erase alpha'ss and cc's, which are not needed in the contraction
!     (just to make the following output clean)
         do j=1,3
            if(j.gt.nmb_cc(i)) then 
               alpha (j,i)=0
               cc(j,i)=0
            endif
         enddo

!     and print everything
         write(6,100) i,center(i),poly(1,i),poly(2,i),poly(3,i),  
     &         nmb_cc(i),(alpha(j,i),j=1,3),(cc(j,i),j=1,3)
      enddo
      
      


      end
