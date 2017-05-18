!-- AT ----------------------------------------------------------------C
!                                                                      C
! THIS LITTLE SAMPLE APPLICATION CALCULATES THE MP2 CORRECTION TO      C
! THE GROUND STATE ENERGY OF A CLOSED-SHELL SYSTEM.                    C
!                                                                      C
! IT SERVES THREE PURPOSES:                                            C
! - TO DEMONSTRATE THE USE OF THE PHIS LIBRARY,                        C
! - TO QUICKLY CHECK NEW BACKENDS OF THIS LIBRARY, AND                 C 
! - TO GIVE BEGINNERS IN THE FIELD OF QUANTUM CHEMISTRY A FEELING OF   C
!   HOW A TYPICAL CALCULATION WORKS.                                   C
!                                                                      C
! THE SAME THING IN C CAN BE FOUND IN SAMPLE.C.                        C
!                                                                      C
!---------------------------------------------------------------- AT --C

      PROGRAM SAMPLE

      IMPLICIT NONE

      INTEGER MOLCAS, GUK, BACKEND
      PARAMETER (MOLCAS=0, GUK=1)

      INTEGER MAXSYM, MAXBAS
      PARAMETER (MAXSYM=8, MAXBAS=1024) 

      INTEGER CAP, LEN

      INTEGER NSYM, NBAS, NATOMS

      INTEGER I, J
      INTEGER A, B
      INTEGER P

      INTEGER IORB, JORB
      INTEGER AORB, BORB
 
      INTEGER ISYM, JSYM
      INTEGER ASYM, BSYM
      INTEGER PSYM
      INTEGER IXJ

      REAL*8 IAJB, IBJA, EIJAB, MP2

      REAL*8 EPSI(MAXBAS), EHF 
      REAL*8 OCC(MAXBAS)
      INTEGER SYM(MAXBAS)

      INTEGER NOCC(8),NVIR(8)
      INTEGER LOCC(MAXSYM,MAXBAS), LVIR(MAXSYM,MAXBAS)

      INTEGER MAJOR,MINOR,PATCH

!...EXTERNAL FUNCTIONS...

      REAL*8 VPQRS
      EXTERNAL VPQRS

      INTEGER PHIS_INIT
      EXTERNAL PHIS_INIT

      INTEGER MT(8,8)
      COMMON /MULTAB/ MT

!...FORMATS...

 100  FORMAT (A4,T8,A3,T16,A3,T24,A)
 110  FORMAT (I4,T8,I3,T16,F3.1,T24,F12.6)
 120  FORMAT (A,F20.12)

!----------------------------------------------------------------------C

! TO BE PORTABLE WE HAVE TO GET THE BACKEND INTERACTIVELY *SIGH*
      WRITE(6,*) "CHOOSE A BACKEND:"
      WRITE(6,*) "MOLCAS - ", MOLCAS
      WRITE(6,*) "GUK - ", GUK
      READ(5,*) BACKEND

      CAP = PHIS_INIT(BACKEND)
      
      IF (.NOT.BTEST(CAP,0).OR..NOT.BTEST(CAP,1).OR..NOT.BTEST(CAP,2)) STOP "BACKEND INCOMPLETE"

      CALL PHIS_GET_INFO( NSYM, NBAS, NATOMS )

      WRITE(6,'(A,I1)') 'NO. OF IRREPS: ', NSYM
      WRITE(6,'(A,I4,/)') 'NO. OF BASIS FUNCTIONS: ', NBAS

      IF (NBAS.GT.MAXBAS) STOP "TOO MANY BASIS FUNCTIONS."

!----------------------------------------------------------------------C

      LEN = MAXBAS
      CALL PHIS_GET_EPSI( EHF, EPSI, LEN )
      CALL PHIS_GET_SYM( SYM, LEN )
      CALL PHIS_GET_OCC( OCC, LEN )

      IF (LEN.LT.0) THEN
! WE CHECKED FOR MAXBAS ABOVE SO THIS SHOULD NEVER HAPPEN
         WRITE(6,*) 'IMPOSSIBLE ERROR: CONTACT THE AUTHOR.'
         STOP
      ENDIF

      WRITE(6,100) 'NO.', 'SYM', 'OCC', 'ORBITAL ENERGY'
      WRITE(6,'(50("-"))')
      DO P = 1,NBAS
         WRITE(6,110) P, SYM(P), OCC(P), EPSI(P)
      END DO

      WRITE(6,'(/,"ELECTRONIC HF-SCF ENERGY: ",F12.6,/)') EHF

!----------------------------------------------------------------------C
      
      DO P = 1,NBAS
         PSYM = SYM(P)
         IF (OCC(P).EQ.2) THEN
            NOCC(PSYM) = NOCC(PSYM) + 1
            LOCC(PSYM,NOCC(PSYM)) = P
         ELSE IF (OCC(P).EQ.0) THEN
            NVIR(PSYM) = NVIR(PSYM) + 1
            LVIR(PSYM,NVIR(PSYM)) = P
         ELSE
            WRITE(6, &
                 '("ORBITAL ",I4," HAD OCCUPATION ",F12.6)') P, OCC(P)
            STOP "THIS IS ONLY FOR CLOSED-SHELL SYSTEMS."
         END IF
      END DO

      DO PSYM = 1,NSYM
         WRITE(6,'("IRREP = ",I1)') PSYM
         WRITE(6,'("OCC: ",10I4)') (LOCC(PSYM,I),I=1,NOCC(PSYM))
         WRITE(6,'("VIR: ",10I4)') (LVIR(PSYM,A),A=1,NVIR(PSYM))
         WRITE(6,*)
      END DO

!----------------------------------------------------------------------C

      CALL PHIS_LOAD_VPQRS

      MP2 = 0
      DO ISYM = 1,NSYM
         DO JSYM = 1,NSYM
            IXJ = MT(ISYM,JSYM)
            DO ASYM = 1,NSYM
               BSYM = MT(IXJ,ASYM)

               DO I = 1,NOCC(ISYM)
                  IORB = LOCC(ISYM,I)
                  DO J = 1,NOCC(JSYM)
                     JORB = LOCC(JSYM,J)
                     DO A = 1,NVIR(ASYM)
                        AORB = LVIR(ASYM,A)
                        DO B = 1,NVIR(BSYM)
                           BORB = LVIR(BSYM,B)
                           IAJB = VPQRS(IORB,AORB,JORB,BORB)
                           IBJA = VPQRS(IORB,BORB,JORB,AORB)
                           EIJAB = EPSI(IORB) + EPSI(JORB) &
                                 - EPSI(AORB) - EPSI(BORB)
                           MP2 = MP2 + (2*IAJB*IAJB - IAJB*IBJA) / EIJAB
                        END DO
                     END DO
                  END DO
               END DO

            END DO
         END DO
      END DO

      WRITE(6,120) 'MP2 CORRECTION:', MP2
      END
