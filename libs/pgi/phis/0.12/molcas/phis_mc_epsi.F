      Subroutine phis_mc_epsi( Ehf, e, len )
      
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...constants...
      Integer LuRun
      Character*8 Lbl
      Parameter ( Lbl = 'OrbE')

C...arguments...
      Real*8 Ehf
      Integer*4 len
      Real*8 e(len)

C...local scalars...
      Integer iRc, iOpt, tmp
      Integer i, j, p
      Integer idx, cnt

C...local arrays...
      Real*8 a(nBasSum)
C      Integer Dummy(MxSym)

*---------------------------------------------------------------------*

#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_epsi_"
#endif
C---  this helps debugging if disc I/O fails
      Call qEnter( 'EPSI' )

C---  space check
      If( nActSum.Gt.len ) Then
         len = -nActSum
         Goto 99
      End If

      Ehf = Escf - PotNuc

      call NameRun(RunFileF)
      LuRun=42
      iRc=1
      iOpt=0
      Call OpnRun(iRc,LuRun, iOpt)
      call Get_dArray(Lbl,a,nBasSum)                 

      Do i = 1,len
         e(i) = 0.0d0
      End Do

C---  Note: if basis functions are deleted due to linear dependencies only 
C---  nOrb orbital energies are calculated and stored. Since MapToActive() 
C---  refers to myBas and we want to end up with nAct orbitals, we need quite 
C---  some counters to keep track of the various indexing schemes:
C---  tmp: basis function, p: orbital, cnt: active orbital,
C---  idx: orbital index as seen by the application, 
      cnt = 0
      tmp = 0
      p = 0
      Do i = 1,nIrrep
         Do j = 1,myBas(i)
            tmp = tmp + 1
            If (j.le.nOrb(i)) Then
               p = p + 1
               idx = MapToExternal(tmp)
               If (idx.ne.0) Then
                  cnt = cnt + 1
                  e(idx) = a(p)
               End If
            End If
         End Do
      End Do

      If (cnt.Ne.nActSum ) Then
         write(0,*) cnt, nActSum
         STOP 'mc_epsi: Inconsistency detected -  contact the author.'
      End if
      len = nActSum

 99   Continue
      Call qExit( 'EPSI' )
      Return

      End subroutine phis_mc_epsi
