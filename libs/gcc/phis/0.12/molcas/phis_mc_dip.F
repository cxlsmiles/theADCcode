      Subroutine phis_mc_dip(X, Y, Z, dim)

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...parameters...
      Integer LuOne
      Character*8 SecName
      Parameter (LuOne = 43, SecName = 'Mltpl  1')

C...arguments...
      Integer*4 dim
      Real*8 x(dim,dim), y(dim,dim), z(dim,dim)

C...local scalars...
      Integer iRc, iOpt


#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_dip_"
#endif

      iRc    = 0

      Call qEnter('mc_dip')

      If (dim.lt.nActSum) Then
         dim = -nActSum
         Goto 999
      End If

      iOpt = 0
      Call OpnOne(iRc, iOpt, OneIntF, LuOne)
      If (iRc.Ne.0) Stop 'Error opening ONEINT file.'

C     x component
      Call mc_proc_one(1, X, nActSum, nBasSum, SecName)

C     y component
      Call mc_proc_one(2, Y, nActSum, nBasSum, SecName)

C     z component
      Call mc_proc_one(3, Z, nActSum, nBasSum, SecName)

      Call ClsOne(iRc, iOpt)

      dim = nActSum

 999  Continue
      Call qExit('mc_dip')
      Return

      End subroutine phis_mc_dip

*======================================================================*
