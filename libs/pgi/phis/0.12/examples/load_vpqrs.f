      PROGRAM GET_INTEGRAL
      IMPLICIT NONE

      INTEGER MOLCAS, GUK, BACKEND,MAXSYM, MAXBAS,ND,NN
      PARAMETER (MOLCAS=131072, GUK=1,MAXBAS=1024)
      PARAMETER (ND=MAXBAS**2)

      INTEGER LEN,CAP
      INTEGER NSYM, NBAS, NATOMS
      INTEGER I,J,K,L,IXJ,IJ,IPOL
      REAL*8 EHF,EPSILON(MAXBAS),OCC(MAXBAS)
      REAL*8 XX(ND),YY(ND),ZZ(ND),DLENGTH
      REAL*8 PX(ND),PY(ND),PZ(ND),DVELOCITY
      REAL*8 INTEGRAL,VPQRS
      INTEGER SYM(MAXBAS),MT(8,8)
      COMMON /MULTAB/ MT

      CHARACTER*3 IRREPS(8)
     &     /'Ag ','B3u','B2u','B1g','B1u','B2g','B3g','Au '/

C...EXTERNAL FUNCTIONS...
      INTEGER PHIS_INIT
      EXTERNAL PHIS_INIT
      EXTERNAL VPQRS

      BACKEND=MOLCAS

C     First step: initializing the PHIS interface and check if the
C     backend provides the functions we need
      CAP = PHIS_INIT(BACKEND)
      PRINT*, 'CAP=',CAP

      IF (.NOT.BTEST(CAP,0) .OR.
     &    .NOT.BTEST(CAP,1) .OR.
     &    .NOT.BTEST(CAP,2) ) STOP "BACKEND INCOMPLETE" ! OR.
C     &    .NOT.BTEST(CAP,11).OR.
C     &    .NOT.BTEST(CAP,12)) STOP "BACKEND INCOMPLETE"

C     Second step: read the 'info' which tells us the nmb of syms, the
C     nmb of basis functions and the nmb of atoms
      CALL PHIS_GET_INFO(NSYM, NBAS, NATOMS)
      IF (NBAS.GT.MAXBAS) STOP "TOO MANY BASIS FUNCTIONS."

C     Reading orbital energies symmetries and occupation numbers and
C     print a table of them to stdout
      LEN = MAXBAS
      CALL PHIS_GET_EPSI( EHF, EPSILON, LEN )
      CALL PHIS_GET_SYM( SYM, LEN )
      CALL PHIS_GET_OCC( OCC, LEN )

      DO I=1,LEN
         PRINT*, I,SYM(I),OCC(I),EPSILON(I)
      END DO
      PRINT*, ' '

C     read all the dipole length moments in the memory
C      NN=MAXBAS
C      CALL PHIS_GET_DIP(XX, YY, ZZ, NN)

C     read all the dipole length moments in the memory
C      NN=MAXBAS
C      CALL PHIS_GET_VEL(PX, PY, PZ, NN)

C print the moments
C       PRINT*,'LENGTH'
C       PRINT*, 'NN = ',NN
C      DO I=1,LEN
C         DO J=1,LEN
C            IJ=(I-1)*NN+J
C            PRINT2002, I,J
C            PRINT2003, IRREPS(SYM(I)),IRREPS(SYM(J))
C            PRINT2004, XX(IJ),YY(IJ),ZZ(IJ)
C            PRINT*, ' '
C         END DO
C      END DO
C      PRINT*,'VELOCITY'
C      PRINT*, 'NN = ',NN
C      DO I=1,LEN
C         DO J=1,LEN
C            IJ=(I-1)*NN+J
C            PRINT2002, I,J
C            PRINT2003, IRREPS(SYM(I)),IRREPS(SYM(J))
C            PRINT2004, PX(IJ),PY(IJ),PZ(IJ)
C            PRINT*, ' '
C         END DO
C      END DO 
C 2002 FORMAT(5X,I4,5X,I4)
C 2003 FORMAT(6X,A3,6X,A3)
C 2004 FORMAT(2X,E12.6,2X,E12.6,2X,E12.6)
C      GO TO 20

C     read all the integrals in the memory
      CALL PHIS_LOAD_VPQRS

C     read the orbital indices 
      PRINT*, 'ENTER THE 4 ORBITAL INDICES'

 10   CONTINUE

      READ(5,*) I,J,K,L

      IF (I.EQ.999) GO TO 20 

C      IF (IPOL.EQ.1) THEN
C         DLENGTH=XX((I-1)*NN+J)
C         DVELOCITY=PX((I-1)*NN+J)
C      ELSE IF  (IPOL.EQ.2) THEN
C         DLENGTH=YY((I-1)*NN+J)
C         DVELOCITY=PY((I-1)*NN+J)
C      ELSE IF (IPOL.EQ.3) THEN
C         DLENGTH=ZZ((I-1)*NN+J)
C         DVELOCITY=PZ((I-1)*NN+J)
C      END IF

C      PRINT*, 'THE LENGTH INTEGRAL =',DLENGTH
C      PRINT*, 'THE VELOCITY INTEGRAL =',DVELOCITY

C    CHECK THE SYMMETRY
      PRINT*, 'THE ORBITAL SYMMETRIES: I J K L'
      PRINT*, SYM(I),SYM(J),SYM(K),SYM(L)
      IXJ = MT(SYM(I),SYM(J))
      IF (SYM(L).NE.MT(IXJ,SYM(K))) THEN
         PRINT*, 'THE INTEGRAL IS ZERO'
      END IF

C     EXTRAXT THE INTEGRAL
      INTEGRAL=VPQRS(I,J,K,L)
      PRINT*, 'THE INTEGRAL =',INTEGRAL
      GO TO 10 

 20   CONTINUE

      STOP
      END
