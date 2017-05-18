      Subroutine phis_mc_nextint(i,j,k,l,Vijkl)

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

      Integer ntibuf
      Parameter (ntibuf = 9600 )
C     If any integral in your system has this value, you loose.
      Real*8 UNSET
      Parameter (UNSET = 12345.6789d0)

C====================== Variable Declarations =========================!

C...Arguments...
      Integer*4 i,j,k,l
      Real*8  Vijkl

C...Scalars...

C nsp is special as it is the outermost loop counter.
C It is initialized in m_setup.
      Integer nsp, nsq, nsr, nss
      Save nsq, nsr, nss
      Common /NSP/ nsp

      Integer nt, nu, nv, nx
      Save nt, nu, nv, nx

      Integer iSec, iOut, OldBlock
      Save iSec, iOut, OldBlock

      Integer nop,noq,nor,nos
      Integer numin,numax
      Integer nxm,ntm,nssm
      Integer nspq,nspqr, nActP
      Integer x

C...Arrays...
      Real*8 tibuf(ntibuf)
      Integer ToC(106),LuTra
      Integer SumBas(0:MxSym), SumAct(0:MxSym)
      Save tibuf, ToC, SumBas, SumAct
      Integer*8 KLI 
      character(len=31) border
      Integer*4 MUL(8,8)
      Common /MULTAB/ MUL

C====================== Executable Statements =========================C


!#ifdef DEBUG
!!VVP:Tracing 
!      print *,"phis_mc_nextint_"
!#endif

      border=repeat("-",31)

 100  Continue

C---  terminate and reset to initial state ---!
      If( nsp.Gt.nIrrep ) Then    
         Call DaClos(20)
         nsp = 0
         l = 0
         Return
      End If

C---  initialize ---!
      If( nsp.Eq.0 ) Then      
#ifdef DEBUG
            Write(6,'(A,/)') 'MapToExternal in phis_mc_nextint:'
            Write(6,'(I5,"->",I5)') (x, MapToExternal(x),x=1,nBasSum)
#endif
        LuTra=20
         Call DANAME(LuTra,TraIntF)

C     get table of contents from TraInt
         iSec = 0
         Call iDAFILE(LuTra,2,ToC,106,iSec)
#ifdef DEBUG
        WRITE(6,'(A)') "" 
        WRITE(6,'(A)') border 
        WRITE(6,'(A)') "*Table of Contents from TraInt*"
        WRITE(6,'(A)') border
        WRITE(6,'(8X,A,17X,A)') "I","IDISK"
        WRITE(6,'(A)') border
        i=1
        do while(ToC(i)/=0)
        WRITE(6,'(5X,I4,10X,I12)') i,ToC(i)
        i=i+1
        enddo
        WRITE(6,'(A)') border
        WRITE(6,'(A)') "" 
#endif
         SumBas = 0
         SumAct = 0

         Do x = 1,nIrrep
            SumBas(x) = SumBas(x-1) + myBas(x)
            SumAct(x) = SumAct(x-1) + nAct(x)
         End Do
CVVP: I am trying to understand situation
         
         nsp = 1 
         nsq = 1
         nsr = 1
         nss = 1

         nv = 1 
         nx = 1
         nt = 1
         nu = 1

         iOut = 0
         OldBlock = 0

      End If

      Vijkl = UNSET

C---  Calculating these at each iteration causes some overhead,
C---  but the code is much more readable.
      nop = myBas(nsp) 
      NOQ=myBas(NSQ)
      NOR=myBas(NSR)
      NOS=myBas(NSS)

      NSPQ=MUL(NSP,NSQ)
      NSPQR=MUL(NSPQ,NSR)
  
      NSSM=NSR
      If(NSR.Eq.NSP) NSSM=NSQ
     
      NUMAX=NOQ
      If(NSP.Eq.NSQ)NUMAX=NT

      NTM=1
      If(NSP.Eq.NSR) NTM=NV

      NXM=NOS
      If(NSR.Eq.NSS) NXM=NV

      If(NSS.Ne.NSPQR) Goto 2000


      nActP = nAct(NSP)*nAct(NSQ)*nAct(NSR)*nAct(NSS)
      If( nActP.Eq.0 ) Goto 2000
     
      If( 1000*nsp+100*nsq+10*nsr+nss.Ne.OldBlock ) Then

         OldBlock = 1000*nsp+100*nsq+10*nsr+nss

#ifdef DEBUG
      write(6,'(a20,4i4)') 'Symmetry block: ',nsr, nss, nsp, nsq
#endif
C         Call dDAFILE(20,2,TIBUF,RtoI*NTIBUF,iSec)
         Call dDAFILE(20,2,TIBUF,RtoI*NTIBUF,iSec)
#ifdef DEBUG
         write(6,*) 'IOUT=', iOut
#endif
         iOut = 0
#ifdef DEBUG
         write(6,*) 'Reading from Section ', iSec
#endif

      End If

     
C---  This is where the work is done. ---!

      K = MapToExternal(SumBas(NSR-1)+NV)
      L = MapToExternal(SumBas(NSS-1)+NX)
      I = MapToExternal(SumBas(NSP-1)+NT)
      J = MapToExternal(SumBas(NSQ-1)+NU)
      KLI=K*L*I

      If(KLI*J.Ne.0) Then
        iOut  = iOut  + 1
        If(IOUT.Gt.NTIBUF) Then
#ifdef DEBUG
          write(6,*) 'Reading from Section ', iSec
#endif
C            Call dDAFILE(20,2,TIBUF,RtoI*NTIBUF,iSec)
          Call dDAFILE(20,2,TIBUF,NTIBUF,iSec)
          IOUT=1
          End If
     
         Vijkl = TIBUF(iOut)

      End If

C---  Now step the counters for the orbitals ---!
      nu = nu + 1
      If( nu.Le.numax ) Goto 9999

      nt = nt + 1
C     numin depends on nt,nv and nx. We have to reset numin (and nu) each time
C     one of these changes.
      NUMIN=1
      If(NSP.Eq.NSR.And.NT.Eq.NV)NUMIN=NX
      nu = NUMIN
      If( nt.Le.NOP ) Goto 9999

      nx = nx + 1
      nt = NTM
      NUMIN=1
      If(NSP.Eq.NSR.And.NT.Eq.NV)NUMIN=NX
      nu = NUMIN
      If( nx.Le.NXM ) Goto 9999

      nv = nv + 1
      nx = 1
      NTM=1
      If(NSP.Eq.NSR) NTM=NV
      nt = NTM
      nu = 1 
C     no need to recalculate numin explicitly as nx=1 anyway.
      If( nv.Le.NOR ) Goto 9999

C---  Step the counters over symmetry blocks ---!
 2000 Continue

      nss = nss + 1
      nv = 1
      nt = 1 
C     no need to recalculate NTM explicitly as nv=1 anyway.
      If( nss.Le.NSSM ) Goto 9999

      nsr = nsr + 1
      nss = 1
      If( nsr.Le.NSP ) Goto 9999

      nsq = nsq + 1
      nsr = 1
      If( nsq.Le.nsp ) Goto 9999

      nsp = nsp + 1
      nsq = 1

 9999 Continue

      If( Vijkl.Eq.UNSET ) Goto 100

      Return

      End subroutine phis_mc_nextint

!Just rewind

      Subroutine phis_mc_nextint_close ()

      Integer nsp
      Common /NSP/ nsp

      logical lod

      inquire (unit=20, OPENED=lod)
      if (lod) call daclos(20)
      if (nsp .ne. 0) nsp = 0
      return

      End subroutine phis_mc_nextint_close



