      Subroutine phis_mc_overlap(S, dim)

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...parameters...
      Integer LuOne
      Character*8 SecName
      Parameter (LuOne = 43,  SecName = 'Mltpl  0')

C...arguments...
      Integer*4 dim
      Real*8 S(dim,dim)

C...local scalars...
      Integer iRc, iOpt

#ifdef DEBUG
!VVP:Tracing 
      print *,"phis_mc_overlap_"
#endif
      iRc    = 0

      Call qEnter('overlap')

      If (dim.lt.nActSum) Then
         dim = -nActSum
         Goto 999
      End If

      iOpt = 0
      Call OpnOne(iRc, iOpt, OneIntF, LuOne)
      If (iRc.Ne.0) Stop 'Error opening OneInt file.'

C     Reads and transforms the overlap integrals.
      Call mc_proc_one(1, S, nActSum, nBasSum, SecName)

      Call ClsOne(iRc, iOpt)

      dim = nActSum

 999  Continue
      Call qExit('overlap')
      Return

      End subroutine phis_mc_overlap

*======================================================================*
