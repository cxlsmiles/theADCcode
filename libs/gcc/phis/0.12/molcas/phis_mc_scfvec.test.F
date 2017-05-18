      Subroutine phis_mc_scfvec( C, n, len )
C
C some notes:
C
C nBasSum is the total number of basis functions, 
C         i.e. the length of an SCF vector
C nActSum is the total number of SCF vectors, which may be smaller than
C         nBasSum due to freezing or deleting orbitals
C nAO is the number of atomic orbitals 
C     (Should be equal to nBasSum? Beware nBas/nOrb confusion!)
C
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...arguments...
      Integer*4 len
      Integer*4 n
      Real*8 C(len,n)

C...local variables...
      Integer*4 nAO, ncc_max
      Real*8 Cso(len,n)

*----------------------------------------------------------------------*

#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_scfvec_"
#endif
      Call qEnter( 'scfvec' )

C---  query for dimensions of the AO basis
      nAO = 0
      ncc_max = 0
      Call phis_mc_ao(0, 0, 0, 0, 0, nAO, ncc_max)
      nAO = -nAO
      ncc_max = -ncc_max
     
C     check boundaries
      If (n.lt.nActSum) n = -nActSum
      If (len.lt.nAO) len = -nAO
      If (n.lt.0 .or. len.lt.0) GOTO 999

      If (len.Ne.nAO) Then
*     This is because I don't pass the leading dimension of C to
*     generate_SOtoAO below.
         Write(0,*) 'FIXME: only works when len = nAO.'
         Stop 'mc_scfec'
      End If

      If (nBasSum.Ne.nAO) Then
         Write(0,*) "Panic: number of basis function is not equal"
         Write(0,*) "to the number of GTOs."
         Write(0,*) "How can that be??"
         Call mc_stop
      End If

      Call mc_read_scfvec(Cso, n, len)

      Call generate_SOtoAO(C, Cso, nActSum, nBasSum, nAO, ncc_max)

      n = nActSum
      len = nAO

 999  Continue
      Call qExit( 'scfvec' )
      Return

      End subroutine phis_mc_scfvec

*======================================================================*

      Subroutine generate_SOtoAO(Dst, Src, nMO, nSO, nAO, ncc_max)
C
C     Generates the mapping from symmetry adapted orbitals (SO) to atomic 
C     orbitals (AO). The first two parameters are passed to allow for 
C     some form of dynamic allocation.
C
C     Short summary, since this one is a bit tricky:
C     * outer loop picks the next SO.
C       Then we know the angular momentum, the number of constituent AOs 
C       and their phases, and the number of SOs in the same shell.
C     * inner loop over the shell size.
C       Locate the next AO in this shell (for a fixed center). Also
C       increment SO counter.
C     * innermost loop over equivalent centers.
C       Build the list of AOs that make up the SO under consideration.
C
C Note: A shell in this context means all functions (on a unique center)
C       with the same angular momentum. This is different from the meaning
C       in quantum chemistry textbooks.
C
C Some important variables:
C SOcenters(i,j), SOphases(i,j): center in the i'th symmetry class of
C                                SO number j, and its phase
C AOcenter(i): center of the i'th AO
C FirstAO(i): offset of the i'th center into the list of AOs
C SOtoAO(i,j): number of the i'th constituent AO of the j'th SO,
C              includes the phase information as sign.
C              Note that SOtoAO(0,j) is the total number of AOs that make
C              the j'th SO.
C

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...Arguments...
      Integer*4 nAO, ncc_max
*     due to a dirty hack, only the first nSO elements of each vector in
*     Src are initialized!!
      Real*8 Src(nAO,nMO), Dst(nAO,nMO)

C...(some) local variables...
      Integer iShell,x,y,z,iDeg
      Integer*4 ncc(nAO)
      Integer*4 poly(3,nAO)
      Real*8 cc(ncc_max,nAO)
      Real*8 alpha(ncc_max,nAO)
      Integer*4 AOcenter(nAO)

      Character*8 SOlabel(mxfnc)
      Integer SOCenters(1:8,nSO), SOPhases(1:8,nSO)
      Integer SOtoAO(0:8,nSO)
      Integer FirstAO(nAtom)

C     Degeneracies for cartesian angular momentum functions
      Integer entart(0:4)
      Data entart /1,3,6,10,15/

*----------------------------------------------------------------------*
#ifdef DEBUG
!VVP:Tracing
      print *,"generate_SOtoAO_"
#endif
!Initialization of the SOtoAO matrix
      SOtoAO=0

C---  Get information on the unique centers.

      Call GETMEM('ICENT','ALLO','INTE',ipCent,8*nSO)
      Call GETMEM('IPHASE','ALLO','INTE',ipPhase,8*nSO)
      Call icopy(8*nSO,0,0,iWork(ipPhase),1)
      Call icopy(8*nSO,0,0,iWork(ipCent),1)

      Call SOout(SOlabel,iWork(ipCent),iWork(ipPhase))

       
      Call icopy(8*nSO,0,0,SOphases,1)
      Call icopy(8*nSO,0,0,SOcenters,1)
      Call SOout(SOlabel,SOcenters,SOphases)
#ifdef DEBUG 
      write(6,*) "===============SoLabel==============="
      do i=1, nSO
         write(6,'(I4,2X,A)') i,SoLabel(i)
      end do 

      write(6,*) "===============SoCenters============="
      do i=1,nSO
         write(6,'(8I4)') SoCenters(1:8,i)
      end do 

      write(6,*) "===============SoPhases=============="
      do i=1,nSO
         write(6,'(8I4)') SoPhases(1:8,i)
      end do
      write(6,*) "====================================="
CC      print *,Nbasis
CC      print *,ipVal
CC      print *,"****************"
CC      print *,nCntr
#endif

C---  Get information on the AOs.
      Call phis_mc_ao(poly, ncc, cc, alpha, AOcenter, nAO, ncc_max)

      If (nAO.lt.0 .or. ncc_max.lt.0) Then
         Write(6,*) "Hot damn!", nAO, ncc_max
         Stop
      End If

C---  Since the AOs are stored in a single array, generate an 
C     auxilliary array that hashes the first AO for a given center.
      iTmp = 0
      Do i = 1,nAO
         If (iTmp.ne.AOcenter(i)) Then
            iTmp = AOcenter(i)
            FirstAO(iTmp) = i
#ifdef DEBUG
            Write(6,'(2(1X,A,1X,I4))') "Atom no.", iTmp,": FirstAO", i
#endif
         End If
      End Do
  
C---  the main loop
      iSO = 1
      Do iDumOuter = 1,nSO
#ifdef DEBUG
         Write(6,'(A)')"========================================="
         Write(6,'(1X,A,I4,A,A,A)') "SO no.   ", iSO, " : Type <<", 
     & SOlabel(iSO),">>"
#endif

C---  Get shell information for the next SO.
         Call parse_lbl(SOlabel(iSO), iCnttp, l, x, y, z)
CVVP.
C          iShellSize = nBasis(ipShll(iCnttp)+l)
           iShellSize = nBasis(ipVal(iCnttp)+l)
#ifdef DEBUG
         write(6,'(A,A,I4)') "parsed label:","iCnttp=",iCnttp
         write(6,'(A,A,I4)') "parsed label:","l     =",l
         write(6,'(A,A,I4)') "parsed label:","x     =",x
         write(6,'(A,A,I4)') "parsed label:","y     =",y
         write(6,'(A,A,I4)') "parsed label:","z     =",z
         Write(6,'(A,4X,I7)')   "shell size  :", iShellSize
#endif
C---  Loop over the given shell.
C     Since I locate AOs based on the polynomial, and the polynomial equal for all
C     AOs in a shell, iLastAO keeps track of which part of the AO vector has already
C     been searched.
         iLastAO = 0
      do ideg=1,nCntr(iCnttp)
         Do iDumInner = 1,iShellSize

C---  Locate AO.
            Do iAO = iLastAO+1,nAO
               If (AOcenter(iAO).eq.SOcenters(1,iSO) .and.
     +              poly(1,iAO).eq.x .and.
     +              poly(2,iAO).eq.y .and.
     +              poly(3,iAO).eq.z) Goto 20
            End Do
c     error condition
            Write(6,*) "Cannot locate AO: ", SOcenters(1,iSO),x,y,z
            Stop

 20         Continue
            iLastAO = iAO

C---  Loop over equivalent centers.
C     iAO refers to centers, iFoo to unique centers, i.e. iFoo is the relative 
C     index of a given AO on a certain type of center.
C     Thus, I can add an offset (FirstAO) to iFoo to get the absolute AO index 
C     for each atom. Note that the number of centers is stored in the 0'th element 
C     for each SO.
            iFoo = iAO - FirstAO(SOcenters(1,iSO))
#ifdef DEBUG
            Write(6,*) "dbg: AO/Foo", iAO, iFoo
#endif
            SOtoAO(0,iSO) = 0
            Do k = 1,8
               If (SOcenters(k,iSO).Ne.0) Then
                  SOtoAO(0,iSO) = SOtoAO(0,iSO) + 1
                  SOtoAO(SOtoAO(0,iSO),iSO) = 
     +               (FirstAO(SOcenters(k,iSO))+iFoo) * SOphases(k,iSO)
               End If
            End Do
#ifdef DEBUG
            Write(6,'(A,4X,8I5)') "constituent AOs", 
     &                         (SOtoAO(i,iSO), i=1,SOtoAO(0,iSO))
#endif
            iSO = iSO + 1
            If (iSO.gt.nSO) Goto 100
         end do 
         End Do
      End Do

 100  Continue

C---  FIXME: 
C     I dont like passing parameters down like this. As of now, the SO->AO
C     mapping is only needed here, but it may be practical to provide it
C     as a separate entity later.

#ifdef DEBUG
      write(6,'(A)') "=============The SOtoAO matrix==============="
      do i=1,nSO
         write(6,'(9I4)') SOtoAO(0:8,i)
         write(6,*) ""
      end do
      write(6,'(A)') "============================================="
#endif
      Call mc_transform_scfvec(Dst, Src, SOtoAO, poly, nMO, nSO, nAO)
      Return
      End subroutine generate_SOtoAO

*======================================================================*

      Subroutine parse_lbl(str, iCnttp, l, x, y, z)
C
C     Helper function that parses a Molcas style basis function label
C     for angular momentum quantum numbers.
C     
C     input: str - the label string, something like 'H   2px' or 'CARB4f2+'
C
C     output: 
C     iCnttp - the type of center relative to some Molcas internal order
C     l - the main angular QN, what Molcas people call 'shell'
C     x,y,z - the exponents, compatible with get_ao().
C
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'
*      Implicit None

C...Arguments...
      Character*(8) str
      Integer iCnttp
      Integer l,x,y,z

C...lacal variables...
      Integer i

*----------------------------------------------------------------------*

      x = 0
      y = 0
      z = 0
#ifdef DEBUG
!VVP: Tracing
      print *,"parse_lbl_"
#endif
c     quick hack to determine the type of center
      Do i = 1,Mxdc
         If (str(1:4).eq.LblCnt(i)(1:4)) Goto 10
      End Do

      Write(6,*) "Cannot determine center type: ", str(1:4)
      Stop

 10   Continue
      iCnttp = i

C---  dispatch on QN l
      Read(str(5:8),'(I1)',err=1234) l
c     Molcas uses a 1-based index, but I want the actual QN
      l = l-1
      Goto (1,2) (l+1)

      If (l.lt.2) Then
         Write(6,*) "error dispatching on main QN: ", l
         Stop
      End If

C---  general case: d orbital or higher
      x = -1
      y = l
      Read(str(7:7),'(I1)') z
      If (str(8:8).eq.'-') z = -z
      Return

C---  s orbital: trivial case
 1    Return

C---  p orbital: 3 functions, all cartesian
 2    Continue
      If (ichar(str(7:7))-ichar('y')) 21,22,23

 21   x = 1
      Return
 22   y = 1
      Return
 23   z = 1
      Return

C---  special case: cartesian functions
 1234 Continue
      If (str(5:5).eq.'d') Then
         l = 2
      Else
         l = ichar(str(5:5)) - ichar('f') + 3
      End If
      Read(str(6:),'(3I1)') x, y, z
      Return

      End subroutine parse_lbl

*======================================================================*

      Subroutine mc_unblock_scfvec(D, Buf, n, len)
*
* MOLCAS stores the SCF vectors in symmetry-blocked form. This routine
* unpacks it into a full matrix D taking care of frozen/deleted orbitals
* on the fly.
* n is the number of active orbitals, len is the length of the SCF vectors.
*
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...Arguments...
      Integer*4 n, len
      Real*8 Buf(*), D(len,n)

C...local variables...
C i,j,k : loop vars
C b,a : cumulative basis/active indices
C o,p : offsets

      Integer i,j,k
      Integer b,a
      Integer o,p

*----------------------------------------------------------------------*
#ifdef DEBUG
!VVP:Tracing
      print *,"mc_unblock_scfvec_"
#endif
      Do i = 1,len
         Do j = 1,n
            D(i,j) = 0.0d0
         End Do
      End Do

C---  Loop over (symmetry-blocked) basis functions (i.e. SOs).
C     NOTE: The number SCF vectors in the INPORB file is equal to the
C           number of basis functions. Vectors that have been deleted
C           due to linear dependencies will contain garbage.
      b = 0
      o = 0
      p = 0
      Do i = 1,nIrrep
         Do j = 1,myBas(i)
            b = b + 1
            a = MapToExternal(b)
C---
C     a is the MO index, k is the (blocked) SO index.
C     o+k is the (unblocked) SO index in the full matrix for the given MO.
C     p is the offset into the blocked matrix.

            If (a.ne.0) Then
               Do k = 1,myBas(i)
                  D(o+k,a) = Buf(p+k)
               End Do
            End If

            p = p + myBas(i)
         End Do

         o = o + myBas(i)
      End Do

      Return
      End subroutine mc_unblock_scfvec

*======================================================================*

      Subroutine mc_transform_scfvec(A, M, SOtoAO, poly, nMO, nSO, nAO)
*
* Takes a matrix of SCF coefficients in the SO basis in M and transforms
* them to the AO basis. The result is returned in A.
*
      Implicit None

C...Arguments...
      Integer*4 nMO, nSO, nAO
      Real*8 M(nSO,nMO), A(nAO,nMO)
      Integer SOtoAO(0:8,nSO)
      Integer*4 poly(3,nAO)

C...local scalars...
      Integer*4 iMO, iSO, iAO
      Integer iCnt, counter
      logical passed(nso)
      Real*8 sf

C     MOLCAS has some silly factors in the cartesian functions.
C     It may be that these are sf = Pi_x sqrt(x),
C     where x are the exponentials in the polynomial,
C     but this has been verified only up to f functions, yet.
      Integer x,y,z
      Real*8 sq3, sq5
      Parameter (sq3 = 1.732050807568877d0,
     +           sq5 = 2.23606797749979d0)
      Real*8 silly_factor(0:3)
      data silly_factor /1.d0, 1.d0, sq3, sq5/

C     We may want to switch off the SO->AO transformation for debugging.
      Logical skip_transform
      Parameter (skip_transform = .False.)

#ifdef DEBUG
!VVP:Tracing
      print *,"mc_transform_scfvec_"
#endif
      If (skip_transform) Goto 123
      A=0.0D0
C---
C     a is the MO index, iSO is the SO index.
C     p is the offset into the blocked matrix.

      Do iMO = 1,nMO
!         passed=.false.
         Do iSO = 1,nSO
            Do iCnt = 1,SOtoAO(0,iSO)
#ifdef DEBUG
               write(6,'(A,I4,2X,A,I4)') "iSO",iSO,"iCnt",iCnt
#endif
               iAO = Abs(SOtoAO(iCnt,iSO)) 
!               if(.not.passed(iAO)) then
               x = poly(1,iAO)
               y = poly(2,iAO)
               z = poly(3,iAO)
               If (x.lt.0) Then
C     spherical functions
                  sf = 1.0d0
               Else
C     cartesian functions

            If ((x+y+z).gt.3) Then
               Write(0,*) 'Error:'
               Write(0,*) 'If you want to use cartesian functions'
               Write(0,*) 'with angular momentum higher thon f,'
               Write(0,*) 'you have to take care of the silly factors'
               Write(0,*) 'in the routine unblock_matrix() first!'
               Stop
            End If

                  sf = silly_factor(x) *
     +                 silly_factor(y) *
     +                 silly_factor(z)
               End If
               A(iAO,iMO) = A(iAO,iMO)+
     +              M(iSO,iMO)/DBLE(SOtoAO(0,iSO)) * sf                
!    +              M(iSO,iMO)*DBLE(SOtoAO(iCnt,iSO))/
!     +              DBLE(SOtoAO(0,iSO)) * sf
!!               A(iAO,iMO) = A(iAO,iMO) / DABS((SOtoAO(iCnt,iSO)))
            if (DBLE(SOtoAO(iCnt,iSO))>0.D0) then
            A(iAO,iMO)=A(iAO,iMO)
            else 
            A(iAO,iMO)=-A(iAO,iMO)
            endif
#ifdef DEBUG
            write(6,*)  iao, imo, A(iAO,iMO) 
#endif
            End Do
         End Do
      End Do

#ifdef DEBUG
      do iMO = 1,nMO
         write(6,'(A)')      "=============================="
         write(6,'(A,2X,I4)')"Final MO LCAO for orbital",imo
         do iSO=1,nSO
            write(6,'(F12.8)') A(iSO,iMO) 
         end do 
!!         write(6,'(A)')      "=============================="
      end Do
#endif
 
      

      Return

 123  Continue
      Write(6,*) "Warning: skipping SO->AO transformation."
      If (nSO.Ne.nAO) Then
         Write(6,*) "number of SOs and AOs must be equal in this case."
         Call mc_stop
      End If

      Do iMO = 1,nMO
         Do iSO = 1,nSO
            A(iSO,iMO) = M(iSO,iMO)
            write(6,*)"NEW", iao,imo,A(iAO,iMO) 
         End Do
      End Do

      Return
      End subroutine mc_transform_scfvec

*======================================================================*

      Subroutine mc_read_scfvec(C, n, len)
*
* Reads the SCF vectors and unpacks them into the matrix C.
* Call mc_unblock_scfvec() to take care of MO deletion/rearrangement.
*
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...Paramters...
      Integer LU
      Parameter (LU = 42)

C...Arguments...
      Integer n, len, ndata
      Real*8 C(len,n)
      
C...local scalars...
      Character*80 VecTit
      Integer lOcc,iErr
      Real*8 Occ
C...local arrays...
      Real*8 A(len*len)
#ifdef V64 
      Integer ndata
      Logical found
#endif
#ifdef DEBUG
!VVP:Tracing
      print *,"mc_read_scfvec_"
#endif
* If lOcc is different from zero occupation numbers are read into Occ.
* snarfed from motra/rdcmo.f
      lOcc = 0
      A=0.0d0
!VVP
#ifdef  V6_2
      Call RdVec(FileName,Lu,'C',nIrrep,myBas,myBas,A,
     &Dummy,Dummy,iDummy,VecTit,1,iErr)
#elif   V6_4 
!!!FIXME
      Call qpg_darray('Last orbitals',found,ndata)
      if(found) Call get_darray('Last orbitals',A,ndata)
#endif

      Call mc_unblock_scfvec(C, A, n, len)

      Return
      End subroutine mc_read_scfvec

*======================================================================*
