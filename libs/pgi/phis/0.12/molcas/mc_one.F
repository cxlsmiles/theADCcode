*
* Helper routines that are common to all routines handling 
* one-electron properties.
*
*----------------------------------------------------------------------*

      Subroutine mc_proc_one(iComp, A, nMO, nSO, SecName)
*
* Read the component iComp of the dipole matrix. The file ONEINT must
* already be open.
* Then expands the matrix from packed storage mode to full storage mode
* and finally transforms it from the SO to the MO basis.
* The final result is returned in A which is an array of dimension (dim,dim).
* The argument SecName determines the section in the ONEINT file that
* the integrals should be taken from.
*
* This reason for wrapping the calls to read_and_expand and mc_one_transform
* in this subroutine is that I want a local dynamic array T(nSO,nSO).
*
      Implicit None

C...Arguments...
      Integer iComp, nMO, nSO
      Real*8 A(nMO,nMO)
      Character*8 SecName

C...local scalars...
      Integer iRc, iOpt, iSymlbl
      Integer buflen
      Integer i,j

C...local arrays...
      Real*8 T(nSO,nSO)
#ifdef DEBUG
!VVP:Tracing
      print *,"mc_proc_one_"
#endif
      iSymlbl = 1
      iOpt = 1
      Call RdOne(iRc, iOpt, SecName, iComp, buflen, iSymlbl)
      If (iRc.Ne.0) Then
         Write(0,*) 'mc_proc_one: Error reading buflen.'
         Call mc_stop
      End If

      Call read_and_expand(iComp, T, nSO, buflen, iSymlbl, SecName)
      Call mc_one_transform(A, T, nMO, nSO)

      Return
      End subroutine mc_proc_one

*----------------------------------------------------------------------*

      Subroutine read_and_expand(iComp, D, n, buflen, iSylbl, SecName)
*
* The block-packed dipole matrix is read into a temporary buffer and then
* expanded into the matrix D which is of dimension (n,n).
* buflen is the size needed for the temporary buffer.
* iSyLbl is the irrep of the matrix, needed to determine which blocks are 
* zero by symmetry.
* SecName is the section in the ONEINT file that the integrals are read from.
*
      Implicit Real*8 (a-h,o-z)
      include 'molcas.fh'

C...Arguments...
      Integer n
      Integer buflen, iSylbl
      Real*8 D(n,n)
      Character*8 SecName

C...local scalars...
      Integer iOpt, iRc
      Integer iOff, jOff, len
      Integer i,j,p
      Integer iSym, jSym, ptr

C...local arrays...
      Real*8 TmpBuf(buflen)
#ifdef DEBUG
!VVP:Tracing
      print *,"read_and_expand_"
#endif
      
      iRc = 0
      iOpt = 6
      Call RdOne(iRc, iOpt, SecName, iComp, TmpBuf, iSym)
      If (iRc.Ne.0) Then
         Write(0,*) 'read_and_expand: Error reading matrix.'
         Write(0,*) SecName, iComp, iSym, iRc
         Call mc_stop
      End If

C---  zero out the destination
      Do i = 1,n
         Do j = 1,n
            D(i,j) = 0.0d0
         End Do
      End Do

      ptr = 1
      Do iSym = 1,nIrrep
         Do jSym = 1,iSym

C---        tricky implementation of the multiplication table
            ij = iEor(iSym-1,jSym-1)
C---        iand->kiand?
            If (iAnd(iSyLbl,2**ij).Eq.0 .Or.
     &           myBas(iSym)*myBas(jSym).Eq.0) Goto 100

C---        diagonal blocks are packed
            If (iSym.Eq.jSym) Then
               len = (myBas(iSym)+1)*myBas(iSym)/2
            Else
               len = myBas(iSym)*myBas(jSym)
            End If

C---        Calculate offsets in the destination matrix
            iOff = 0
            Do i = 1,iSym-1
               iOff = iOff + myBas(i)
            End Do

            jOff = 0
            Do j = 1,jSym-1
               jOff = jOff + myBas(j)
            End Do

C---        the actual copy loop is here
            p = 0
            Do j = 1,myBas(jSym)
               If (iSym.eq.jSym) Then
                  iMax = j
               Else
                  iMax = myBas(iSym)
               End If
               Do i = 1,iMax
                  D(iOff+i,jOff+j) = TmpBuf(ptr+p)
                  D(jOff+j,iOff+i) = TmpBuf(ptr+p)
                  p = p + 1
               End Do
            End Do

C---        some error checking
            If (p.Ne.len) Then
               Write(0,*) 'read_and_expand'
               Write(0,*) 'Bug: wrong number of elements.', p, len
               Call mc_stop
            End If
            If (ptr+p-1.Gt.buflen) Then
               Write(0,*) 'read_and_expand'
               Write(0,*) 'Bug: pointer past end of array.', ptr, buflen
               Call mc_stop
            End If

            ptr = ptr + len
 100        Continue
         End Do
      End Do

      Return
      End subroutine read_and_expand

*----------------------------------------------------------------------*

      Subroutine mc_one_transform(M, A, nMO, nSO)
*
* Transforms a one-electron matrix in A from the SO to the MO basis
* and returns the result in M.
* n and nSO are the dimension of the MO and SO basis, respectively.
*  
      Implicit None

C...arguments...
      Integer nMO, nSO
      Real*8 M(nMO,nMO), A(nSO,nSO)

C...local arrays...
      Real*8 T(nSO,nSO), C(nSO,nMO)
#ifdef _DEBUG:
!VVP:Tracing 
      print *,"mc_one_transform_" 
#endif
      Call mc_read_scfvec(C, nMO, nSO)

C     T = A.C
      Call DGeMm('N', 'N', nSO, nMO, nSO, 1.0d0, 
     &           A, nSO, C, nSO, 0.0d0, T, nSO)
C     M = Transpose(C).T
      Call DGeMm('T', 'N', nMO, nMO, nSO, 1.0d0, 
     &           C, nSO, T, nSO, 0.0d0, M, nMO)

      Return
      End subroutine mc_one_transform

C ==================================================================== C
