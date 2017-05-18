      Subroutine mc_rd_runfile
C
C     This routine sets nOrb(8),nAct(8) and nOcc(8) in the common block from
C     COMFILE.
C
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...parameters...
      Integer MxBuf
      Parameter (MxBuf = 16384)
      
C...local scalars...
      Character*16 Lbl
      Integer iRc, iOpt
      Real*8 Etmp
      Integer mSym, i

C...local arrays...
      Integer mBas(MxSym), mOcc(MxSym), mOrb(MxSym),mAct(MxSym)
!      Real*8 Buf(MxBuf)
#ifdef DEBUG
!VVP:Tracing
      print *, "mc_rd_runfile_"
#endif

!Initizialization of arrays
      do i=1,MxSym
      mBas(i)=0
      mOcc(i)=0
      mOrb(i)=0
      mAct(i)=0

      nOrb(i)=0
      nOcc(i)=0
      nAct(i)=0
      
      end do 

      Call qEnter( 'RunFileInfo' )

      call NameRun(RunFileF)
#ifdef DEBUG
      iRc=1
      iOpt=0
      call DumpRun(iRc, iOpt)
#endif
      LuRun=11
      iRc=1
      iOpt=0
      Call OpnRun(iRc,LuRun, iOpt)
      Lbl='Last energy'
      call Get_dScalar(Lbl,Etmp)
      Lbl='nSym'       
      call Get_iScalar(Lbl,mSym)
      Lbl='nOrb'
      call Get_iArray(Lbl,mOrb(1:mSym),mSym)
      Lbl='nBas'       
      call Get_iArray(Lbl,mBas(1:mSym),mSym)
      Lbl='nIsh'       
      call Get_iArray(Lbl,mOcc(1:mSym),mSym)
      Lbl='nAsh'       
      call Get_iArray(Lbl,mAct(1:mSym),mSym)
      call DaClos(LuRun)
      nIrrep=mSym
      Escf = Etmp

      Do i = 1,nIrrep
         nOrb(i) = mOrb(i)
         nOcc(i) = mOcc(i)
         nAct(i) = mAct(i)
      End Do
      Call qExit( 'RunFileInfo' )
      Return

      End subroutine mc_rd_runfile

C----------------------------------------------------------------------C

      Subroutine mc_rd_traone
C
C     This routine sets the arrays nFro(8) and nDel(8) in the 
C     common block from TRAONE.
C
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...parameters...
      Integer LuOne

C...local scalars...
      Integer LblSize
      Parameter (LblSize = (2*MxOrb*4)/ItoB)
      Integer i!, iDum, iDisk
      Integer mSym
      Real*8 ignore

C...local arrays...
C      Integer Toc(1024)
      Integer IoList(64)
      Character*4 BasLbl(2,MxOrb)
      Integer mBas(MxSym)
C      Integer Dummy(MxSym)

#ifdef DEBUG
!VVP: Tracing
      print *, "mc_rd_traone_"
#endif
      Call qEnter( 'TraOne' )

C     There is no high level interface to TraOne as of Molcas 5.0,
C     so we have to use the low level routines.

      LuOne=20  
      Call DANAME(LuOne,TraOneF)

C---  header information
      iDisk = 0
      Call WR_MOTRA_Info(LUONE,2,iDisk,
     &                   IoList,64,
     &                   ignore,
     &                   mSym,
     &                   mBas,nAct,nFro,nDel,8,
     &                   BasLbl,4*2*MxOrb)

CVVP: In any case, we must to get nAct array.
CPlease,see subroutine mc_rd_runfile
C
CVVP:Old interface to 'TRAONE' file. This interface is valid
Cfor molcas v4.*-v5.*
C      Call GSLIST(IoList,8, Toc,64,
C     &     ignore,RtoI, 
C     &     mSym,1, mBas,8, nAct,8, nFro,8, nDel,8, BasLbl,LblSize ) 
C
C      Call DAFILE(LUONE,4,IoList,iDum,iDisk)

      Call DACLOS(LUONE)
C---  sanity checks
      If (nIrrep.Ne.mSym) Stop 'Inconsistent no. of irreps in TRAONE.'
      Do i = 1,nIrrep
         If (myBas(i).Ne.mBas(i)) Then
            Write(6,'(A,3I5)') 'Inconsistency in TRAONE:',
     +           i ,myBas(i), mBas(i)
            Stop
         End If

         If ((myBas(i)-nDel(i)).Gt.nOrb(i)) Then
            Write(6,'(A,4I5)') 
     +      'You must delete all linear dependent orbitals in MOTRA.',
     +      i, myBas(i), nDel(i), nOrb(i)
            Stop
         End If
      End Do

      Call qExit( 'TraOne' )
      Return

      End subroutine mc_rd_traone

C----------------------------------------------------------------------C
C Note: As I am linking the MA library anyway, I might just as well
C       allocate WorkSpc() locally.
C
C This routine is called from the C init hook to set up all kinds of
C (FORTRAN) internal data. At the moment this is basically just the
C functions for mapping the orbitals.
C
      Subroutine mc_setup(Flags, WorkSpc, wLen)

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...arguments...
      Integer*4 Flags
      Integer*4 wLen
      Real*8 WorkSpc(wLen)

C...local variables...
      Integer Irrep, i
      Integer ipInf
      Integer iBas_Lab(mxdc)
      Character*4 AtomLabel(mxdc)

C This is for nextInt
      Integer nsp
      Common /NSP/ nsp

#ifdef DEBUG
!VVP:Tracing
      print *,"mc_setup_"
#endif
      nsp = 0

C--- Initialize seward data structures, see MOLCAS/src/scf/interf_scf.f

      call NameRun(RunFileF)
      Call Inter1(AtomLabel,iBas_Lab,Coor,Znuc,nAtom,ipInf)
C--- nBas is defined to start from 0 in info.inc, so I use my own array
C      call DaClos(LuRun)
      Do i = 1,nIrrep
         myBas(i) = nBas(i-1)
#ifdef DEBUG 
      write(6,'(2(2X,A,I4))') "Irrep Nr. ",i, "nBas:",myBas(i)
#endif     
      End Do

      Call mc_rd_runfile()
      Call mc_rd_traone()

#ifdef DEBUG
      Write(6,*) 'DEBUG'
      Write(6,'(A,8I5)') 'nBas', (myBas(i), i=1,nIrrep)
      Write(6,'(A,8I5)') 'nOrb', (nOrb(i), i=1,nIrrep)
      Write(6,'(A,8I5)') 'nDel', (nDel(i), i=1,nIrrep)
      Write(6,'(A,8I5)') 'nFro', (nFro(i), i=1,nIrrep)
      Write(6,'(A,8I5)') 'nAct', (nAct(i), i=1,nIrrep)
      Write(6,*) 'END DEBUG'
#endif
C      nBasSum = sum(myBas)
C      nActSum = sum(nAct)


      Do Irrep = 1,nIrrep
         nBasSum = nBasSum + myBas(Irrep)
         nActSum = nActSum + nAct(Irrep)
      End Do

      If (nBasSum.Gt.MxBas) STOP 'Out of memory in m_setup.'
      Call mc_maps(nBasSum, Flags)

      Return
      End subroutine mc_setup

C----------------------------------------------------------------------C

      Subroutine mc_maps( len, Flags )
C
C     This routine generates the array MapToExternal() that maps the 
C     internal orbital indices to the ones seen by the application.
C     Indices can differ in two ways: orbitals may have been deleted
C     during the MO transformation and/or the application may have 
C     requested an energy-first ordering of orbitals.
C
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...arguments...
      Integer*4 len
      Integer*4 Flags

C...local variables...
      Real*8 dum
      Integer cnt, act, tmp
      Integer iErr,Two,Zero
      parameter (Two=2,Zero=0)
      Integer MapToActive(len), MapUnblock(0:len)
      Real*8 TmpSpc(len)

C---  1st step: orbital -> active orbital mapping
C     Inactive functions are marked as zero.

#ifdef DEBUG
!VVP:Tracing
      print *,"mc_maps_"
#endif
      Do i = 1,len
         MapToActive(i) = 0
      End Do

      cnt = 0
      act = 0
      Do Irrep = 1,nIrrep
         Do i = 1,myBas(Irrep)
            cnt = cnt + 1
            If (i.Gt.nFro(Irrep) .and. 
     &           i.Le.(myBas(Irrep)-nDel(Irrep))) Then
               act = act + 1
               MapToActive(cnt) = act
            End If
         End Do
      End Do

C---  2nd step: symmetry first -> energy first sorting
CVVP iAnd->KiAnd? Especially for PGI compilers(i.e.,pgf77,pgf90,pgf95)
      If (iAnd(Flags,2).Eq.0) Then
C---  energy-first ordering requested

c        Need to set MapToExternal temporarily, 
c        because mc_phis_epsi() below uses it.
         Do i = 1,len
            MapToExternal(i) = MapToActive(i)
         End Do

         tmp = len
         Call phis_mc_epsi(dum, TmpSpc, tmp)
         If (tmp.Lt.0) Stop "mc_maps: blown away..."
         Call mc_dpsort(TmpSpc, tmp, MapUnblock(1), 1, iErr)
         If (iErr.Ne.0) Stop 'Error in dpsort.'
         Call mc_invert(MapUnblock(1), tmp, TmpSpc)

      Else

         Do i = 1,len
            MapUnblock(i) = i
         End Do

      End If

C---  3rd step: concatenate to get the final map

      MapUnblock(0) = 0
      Do i = 1,len
         MapToExternal(i) = MapUnblock(MapToActive(i))
#ifdef DEBUG
         Write(6,'(A,2X,4I4)') "DBG: ", 
     &        i, MapToActive(i), MapUnblock(i), MapToExternal(i)
#endif
      End Do


      Return
      End subroutine mc_maps

C----------------------------------------------------------------------C

      Subroutine unblock_packed( src, dest )
C
C     Called from get_overlap to rebuild the full matrix
C     FIXME: matrix dimensions should be explicit in the function call
C
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...Arguments...
      Real*8 src(*), dest(nBasSum,nBasSum)

C...local variables...
C i,j,k : loop vars
C b,a : cumulative basis/active indices
C p,o : src/target offsets

      Integer i,j,k
      Integer a
      Integer o,p

#ifdef DEBUG
!VVP:Tracing

      print *,"unblock_packed_"
#endif
      Do i = 1,nBasSum
         Do j = 1,nBasSum
            dest(i,j) = 0.d0
         End Do
      End Do

      a = 0
      o = 0
      p = 0
      Do i = 1,nIrrep
         Do j = 1,myBas(i)
            a = a + 1
            Do k = 1,j
               dest(o+k,a) = src(p+k)
               dest(a,o+k) = dest(o+k,a)
            End Do
            p = p + j
         End Do
         o = o + myBas(i)
      End Do

      Return
      End subroutine unblock_packed

C----------------------------------------------------------------------C

      Subroutine mc_invert(Map, n, Work)
C
C     We want the inverse mapping w.r.t. what dpsort() returns.
C     MapUnblock(i) is defined as the position of the ith orbital in energy 
C     order so we can conveniently write 
C        A(MapUnblock(i)) = value
C     in assignments.
C
      Implicit None

C...Arguments...
      Integer n
      Integer Map(n)
      Integer Work(n)

C...local variables...
      Integer i

#ifdef DEBUG
!VVP:Tracing
      print *,"mc_invert_"
#endif
      Do i = 1,n
         Work(Map(i)) = i
      End Do

      Do i = 1,n
         Map(i) = Work(i)
      End Do
  
      Return
      End subroutine mc_invert

*----------------------------------------------------------------------*

      Function mc_query_vel()
*
* Little helper function that tries to read the velocity integrals and
* returns the error code of that attempt.
* Called by the setup routine when setting the capabilities.
*

      Implicit None

C...      include 'molcas.fh'

C...Parameters...
      Integer LuOne
      Parameter (LuOne = 42)

C...Arguments...
      Integer mc_query_vel

C...local variables...
      Integer iOpt, iRc, iComp, iSym
      Integer dummy
#ifdef DEBUG
!VVP:Tracing
      print *,"mc_query_vel_"
#endif

      iComp = 1
      iSym = 1
      iRc = 0

#ifdef DEBUG
CVVP:Dumping of the one-electron integral file('TRAONE')
      iOpt = 1024
#else
      iOpt=0
#endif
      call NameRun('RUNFILE')
      Call OpnOne(iRc, iOpt, 'ONEINT', LuOne)
      If (iRc.Ne.0) Stop 'Error opening ONEINT file.'

      iOpt = 1
      Call RdOne(iRc, iOpt, 'Velocity', iComp, dummy, iSym)
      mc_query_vel = iRc

      Call ClsOne(iRc, iOpt)

      Return
      End function mc_query_vel

C ==================================================================== C
