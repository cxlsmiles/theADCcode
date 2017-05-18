      Subroutine phis_mc_info( nSym, n, nCenters )

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...arguments...
      Integer*4 nSym, n, nCenters

C...local scalars...
      Integer i

#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_info_"
#endif
      nSym = nIrrep
      n = nActSum
      nCenters = nAtom

      Return

      End subroutine phis_mc_info
