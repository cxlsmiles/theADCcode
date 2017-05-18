      Subroutine phis_mc_occ( o, len )

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...arguments...
      Integer*4 len
      Real*8 o(len)

C...local scalars...
      Integer i, j, off, idx, tmp

*---------------------------------------------------------------------*

#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_occ_"
#endif
C---  space check
      If( nActSum.Gt.len ) Then
         len = -nActSum
         Goto 99
      End If

      Do idx = 1,nActSum
         o(idx) = 0.d0
      End Do

      tmp = 0
      off = 0
      Do i = 1,nIrrep
         Do j = 1,nOcc(i)
            tmp = off + j
            idx = MapToExternal(tmp)
            If (idx.Ne.0) o(idx) = 2.d0
         End Do
         off = off + myBas(i)
      End Do
      len = nActSum

 99   Continue
      Return

      End subroutine phis_mc_occ
