      Subroutine phis_mc_geometry(n, geo, Z_nuc, Enuc)

C FIXME: Z_nuc is not set!

      Implicit Real*8 (a-h,o-z)
      Include 'molcas.fh'

C...arguments...
      Integer n
      Real*8 geo(3,n)
      Real*8 Z_nuc(n)
      Real*8 Enuc

C...local variables...
      Integer i

*======================================================================*

#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_geometry_"
#endif
C---  space check
      If (n.lt.nAtom) Then
         n = -nAtom
         Return
      End If

      Do i = 1,nAtom
         geo(1,i) = Coor(1,i)
         geo(2,i) = Coor(2,i)
         geo(3,i) = Coor(3,i)
      End Do

      Do i = 1,nAtom
         Z_nuc(i) = Znuc(i)
      End Do

      Enuc = PotNuc

      Return
      End subroutine phis_mc_geometry
