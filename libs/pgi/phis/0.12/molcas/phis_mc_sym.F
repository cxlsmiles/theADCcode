      Subroutine phis_mc_sym( s, len ) 

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...arguments...
      Integer*4 len
      Integer*4 s(len)

C...local scalars...
      Integer i, j, idx
      Integer tmp, cnt

*---------------------------------------------------------------------*

#ifdef _DEBUG_
!VVP:Tracing
      print *,"phis_mc_sym_"
#endif
C---  space check
      If( nActSum.Gt.len ) Then
         len = -nActSum
         Goto 99
      End If

      cnt = 0
      tmp = 0
      Do i = 1,nIrrep
         Do j = 1, myBas(i)
            tmp = tmp + 1
            idx = MapToExternal(tmp)
            If (idx.Ne.0) Then
               cnt = cnt + 1 
               s(idx) = i
            End If
         End Do
      End Do

      If( cnt.Ne.nActSum ) Then
         Write(0,*) cnt, nActSum
         STOP 'mc_sym: Inconsistency detected -  contact the author.'
      End If
      len = nActSum

 99   Continue
      Return

      End subroutine phis_mc_sym
