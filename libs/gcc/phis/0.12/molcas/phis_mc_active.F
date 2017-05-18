      Subroutine phis_mc_active( list, n )

      Implicit Real*8 (a-h,o-z)
      include'molcas.fh'

C...Arguments...
      Integer*4 n
      Integer*4 list(n)

C...local variables...
      Integer i

*======================================================================*

#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_active_"
#endif
      If (n.lt.nActSum) Then
         n = -nActSum
         Return
      End If
      
      
      Do i = 1,nBasSum
         If (MapToExternal(i).ne.0) list(MapToExternal(i)) = i 
      End Do

      n = nActSum
      Return
      End subroutine phis_mc_active
