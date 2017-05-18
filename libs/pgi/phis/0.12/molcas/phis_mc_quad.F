      Subroutine phis_mc_quad(XX, XY, XZ, YY, YZ, ZZ, dim)

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...parameters...
      Integer LuOne
      Character*8 SecName
      Parameter (LuOne = 43, SecName = 'Mltpl  2')

C...arguments...
      Integer*4 dim
      Real*8 XX(dim,dim), XY(dim,dim), XZ(dim,dim),
     +       YY(dim,dim), YZ(dim,dim), ZZ(dim,dim)

C...local scalars...
      Integer iRc, iOpt


#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_quad_"
#endif
      iRc    = 0

      Call qEnter('mc_quad')

      If (dim.lt.nActSum) Then
         dim = -nActSum
         Goto 999
      End If

      iOpt = 0
      Call OpnOne(iRc, iOpt, OneIntF, LuOne)
      If (iRc.Ne.0) Stop 'Error opening OneInt file.'

C     get the individual components of the quadrupole moment
      Call mc_proc_one(1, XX, nActSum, nBasSum, SecName)
      Call mc_proc_one(2, XY, nActSum, nBasSum, SecName)
      Call mc_proc_one(3, XZ, nActSum, nBasSum, SecName)
      Call mc_proc_one(4, YY, nActSum, nBasSum, SecName)
      Call mc_proc_one(5, YZ, nActSum, nBasSum, SecName)
      Call mc_proc_one(6, ZZ, nActSum, nBasSum, SecName)

      Call ClsOne(iRc, iOpt)

      dim = nActSum

 999  Continue
      Call qExit('mc_quad')
      Return

      End subroutine phis_mc_quad

*======================================================================*
