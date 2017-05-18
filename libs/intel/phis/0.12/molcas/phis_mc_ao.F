* FIXME: 
* normalize coefficients

      Subroutine phis_mc_ao(polynomial, ncc, cc, alpha, 
     +                      center, nao, ncc_max)

      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C I use these from MOLCAS/include/info.inc
C nCnttp: number of (symmetry unique) centers
C nTot_Shells: the number of shells associated with the i'th basis set
C

C...arguments...
      Integer*4 nao, ncc_max
      Integer*4 ncc(nao)
      Integer*4 polynomial(3,nao)
      Real*8 cc(ncc_max,nao)
      Real*8 alpha(ncc_max,nao)
      Integer*4 center(nao)

C...local variables...
      Integer x,y,z
      Integer iShellSize

*======================================================================*

#ifdef DEBUG
!VVP:Tracing
      print *,"phis_mc_ao_"
#endif

      If (nao.lt.nBasSum) nao = -nBasSum
      
C      call Inter1(AtomLabel,iBas_Lab,Coor,Znuc,nAtom,ipInf)

C---  Loop over the shells once to get the maximum contraction length.
C     No idea how this works, copied it from MOLCAS/src/scf/interf_scf.f
      itmp = 0

      Do iCnttp=1,nCnttp
CVVP: The source is the same [/src/scf/interf_scf.f]
         Do l=0,nVal_Shells(iCnttp)-1
            ishell=l+ipVal(iCnttp)
#ifdef DEBUG
            write(6,*) 'AO dbg:', icnttp, Transf(iShell)
#endif
            Do i = 1,iShell
               if (nExp(i).gt.itmp) itmp = nExp(i)
            End Do
         End Do
      End Do

      If (ncc_max.lt.itmp) ncc_max = -itmp

      If (nao.lt.0 .or. ncc_max.lt.0 ) Return

      itmp = 0
      iAtom = 0
      iAO = 0

C---  Generates a loop over all actual atoms from the information about
C     unique centers and their degeneracy.

      Do iCnttp=1,nCnttp
         Do iCntr = 1,nCntr(iCnttp)
            itmp = itmp + 1
            nDeg = nIrrep/nStab(itmp)
            Do iDeg = 1,nDeg
               iAtom =iAtom + 1

C---  Loops over all functions in all shells for the center under consideration.
CVVP: Old code of TS
C             yyy=nTot_Shells(iCnttp)-1-(nPrj_Shells(iCnttp)+
C     &nSRO_Shells(iCnttp))
C               Do l=0,yyy
C                  ishell=ipVal(iCnttp)+l
               Do l=0,nVal_Shells(iCnttp)-1
                  ishell=l+ipVal(iCnttp)
#ifdef DEBUG
        Write(6,*) 'AO dbg, iShellSize:', nBasis(iShell)
#endif
                  Do iShellSize = 1, nBasis(iShell)
#ifdef DEBUG
        Write(6,*) 'AO dbg:', iAO, iAtom, iShell, nExp(iShell),
     &             (ipExp(iShell)+i, i=0,nExp(iShell)-1),
     &             (Work(ipExp(iShell)+i), i=0,nExp(iShell)-1)
        Write(6,*) 'AO dbg:'
#endif

C---  There is no stored information on the angular momentum, so we generate
C     that here.

C---                 spherical functions
		                                                                                                      
                     If (l.ge.2 .and. Transf(iShell)) Then
                        Do m = -l,l
                           iAO = iAO + 1                     
                           center(iAO) = iAtom

                           polynomial(1,iAO) = -1
                           polynomial(2,iAO) = l
                           polynomial(3,iAO) = m
  
                           ncc(iAO) = nExp(iShell)
                           Do iExp = 1,nExp(iShell)
                              alpha(iExp,iAO) = 
     +                             Work(ipExp(iShell)+iExp-1)
                              cc(iExp,iAO) = Work(ipCff(ishell) + 
     +                             (iShellSize-1)*nExp(ishell)+iExp-1)
                           End Do

                        End Do

C---                 cartesian functions, contains s and p functions as special cases
                     Else
                       Do x = l,0,-1
                           Do y = l-x,0,-1
                              z = l-x-y
                              iAO = iAO + 1                     
                              center(iAO) = iAtom
                              polynomial(1,iAO) = x
                              polynomial(2,iAO) = y
                              polynomial(3,iAO) = z 

C--   FIXME: we are really nesting too deep here.
                              ncc(iAO) = nExp(iShell)
                              Do iExp = 1,nExp(iShell)
                                 alpha(iExp,iAO) = 
     +                             Work(ipExp(iShell)+iExp-1)
                                 cc(iExp,iAO) = Work(ipCff(ishell) + 
     +                             (iShellSize-1)*nExp(ishell)+iExp-1)
                              End Do

                           End Do
                       End Do

                     End If

                  End Do
               End Do
            End Do
         End Do
      End Do

      Return
      End subroutine phis_mc_ao
