      Subroutine phis_mc_oneel(H, T, dim)

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...parameters...
      Integer LuOne
      Character*8 HamName, KinName
      Parameter (LuOne = 43, KinName = 'Kinetic ', HamName = 'OneHam  ')

C...arguments...
      Integer*4 dim
      Real*8 H(dim,dim), T(dim,dim)

C...local scalars...
      Integer iRc, iOpt


#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_oneel_"
#endif
      iRc    = 0

      Call qEnter('mc_oneel')

      If (dim.lt.nActSum) Then
         dim = -nActSum
         Goto 999
      End If

      iOpt = 0
      Call OpnOne(iRc, iOpt,OneIntF, LuOne)
      If (iRc.Ne.0) Stop 'Error opening ONEINT file.'

C     One electron Hamiltonian
      Call mc_proc_one(1, H, nActSum, nBasSum, HamName)

C     Kinetic energy
      Call mc_proc_one(1, T, nActSum, nBasSum, KinName)

      Call ClsOne(iRc, iOpt)

      dim = nActSum

 999  Continue
      Call qExit('mc_oneel')
      Return

      End subroutine phis_mc_oneel

*======================================================================*
