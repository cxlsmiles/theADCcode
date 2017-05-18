
C this file  reads  all input datas for CAP and generates  
C CAP 1-e integrals in AO basis [modified by nico on Sept-22-2010]
! Yasen: the box dimensions are passed now as input parameters
      SUBROUTINE Compute_CAP_MO(XB, YB, ZB, NB,ISPC)     
      IMPLICIT REAL*8 (A-H,O-Z)

      integer *4 basis_n(ISPC)
      integer *4 ang_pow(ISPC,3)
      
      real *8 ord_expon(ISPC),ord_coe(ISPC)
      real *8 R_m(4),R_n(4)
      real *8 R_ms(ISPC,3)
      real *8 BOX(3),RLAMBDA
      real *8 CAP_AO((NB*(NB+1))/2)
C      real *8 CAP_MO((NB*(NB+1))/2)
      real *8 S_AO((NB*(NB+1))/2)
C      real *8 S_MO((NB*(NB+1))/2)
      real *8, allocatable :: AOMO(:,:)
      real *8 core(3),cen(3)
      common /boxsize/ core
      common /symcent/cen
      CHARACTER *8 IA
      EXTERNAL IPOINT
      

      OPEN(unit=14,file='xyz.txt',STATUS='OLD',FORM='FORMATTED')

      READ(14,*) n_atoms

      DO  k=1,n_atoms
         READ(14,*)IA,UW,X,Y,Z
         if(K.eq.1)then
         XP=X
         XM=X
         YP=Y
         YM=Y
         ZP=Z
         ZM=Z
         endif
         if(k.ne.1)then
            if(X.gt.XP)XP=X
            if(X.lt.XM)XM=X
            if(Y.gt.YP)YP=Y
            if(Y.lt.YM)YM=Y
            if(Z.gt.ZP)ZP=Z
            if(Z.lt.ZM)ZM=Z
         endif
      END DO
      close(14)
C      write(*,*)XP,XM
C      write(*,*)YP,YM
C      write(*,*)ZP,ZM
      
      write(*,*)'BOX-X, BOX-Y,BOX-Z'
!      read(*,*)XB,YB,ZB
      
      CEN(1)=(XP+XM)/2.0d0
      CORE(1)=XP-CEN(1)+XB
      CEN(2)=(YP+YM)/2.0d0
      CORE(2)=YP-CEN(2)+YB
      CEN(3)=(ZP+ZM)/2.0d0
      CORE(3)=ZP-CEN(3)+ZB
      write(*,*)CORE(1),CEN(1)
      write(*,*)CORE(2),CEN(2)
      write(*,*)CORE(3),CEN(3)
C      stop
C      write(*,*)'X,X0'
C      read(*,*)CORE(1),CEN(1)
C      write(*,*)'Y,Y0'
C      read(*,*)CORE(2),CEN(2)
C      write(*,*)'Z,Z0'
C     read(*,*)CORE(3),CEN(3)
      ang_pow(:,:) = 0
      CALL INPREAD(N_B,ISPC,ord_expon,ord_coe,R_ms,basis_n,ang_pow)
      write(*,*)N_B,NB
 
C JUST CHECK THAT WE DONT NEED NOW    
C      CALL CONGUANORM(NB,ISPC,ord_expon,ord_coe,R_ms,basis_n,ang_pow)
C    CALL SMATRIX(NB,ISPC,ord_expon,ord_coe,R_ms,basis_n,ang_pow,S_AO) 
 
CNICO 25.02.2011
C       WRITE(*,*)"DO NOT FORGET TO UNCOMMENT AOINTEGRAL"  
      CALL AOINTEGRAL(NB,ISPC,ord_expon,ord_coe,R_ms,basis_n,ang_pow,
     *  CAP_AO)
C       WRITE(*,*)"DO NOT FORGET TO UNCOMMENT AOINTEGRAL"  
!      print *, CAP_AO
!      stop
C JUST CHECK THAT WE DONT NEED NOW    
C      k=0
C      open(unit=21,file='aosmat.dat',STATUS='UNKNOWN',FORM='FORMATTED')
C      open(unit=22,file='aocapmat.dat',STATUS='UNKNOWN',
C     *     FORM='FORMATTED')
C      DO i=1,N_B
C         DO j=1,i
C            k=k+1
C            write(22,*)i,j,CAP_AO(k)
C            write(21,*)i,j,S_AO(k)
C         end do
C      end do
      
      open(unit=8,file='mocoef.txt',STATUS='UNKNOWN',FORM='FORMATTED')
      rewind(8)
C N_B is the number of AOs, NB1 is the number of active MOs
      read (8, *) N_B, NB1
      IF (NB .ne. N_B) THEN
         print *, "Error reading the number of basis functions"
         stop
      ENDIF

      allocate(AOMO(NB1,N_B))
      AOMO(:,:) = 0d0
      
      DO i=1,NB1
         DO j=1,N_B
            read (8, *)ik, kl, AOMO(i,j)
         end do
      end do
      close(8)

      CALL CITR(NB,NB1,AOMO,CAP_AO,'capmat.dat') 
      CALL CITR(NB,NB1,AOMO,S_AO,'mosmat.dat')    
      deallocate(AOMO)

      RETURN
      END
      

C!=============================================================================
      SUBROUTINE INPREAD(N_B,ISPC,ord_expon,ord_coe,R_ms,basis_n,
     *   ang_pow)

      IMPLICIT REAL*8 (A-H,O-Z)
      
      CHARACTER *1 shell_type,BLNK(6)
      CHARACTER *3 shell_p(3),shell_d(6),Shell_f(10)
      CHARACTER *4 Shell_G(15) 
      CHARACTER *3 shell_arr(ISPC)   
      CHARACTER *8 IA
      

      integer *4 angu_pow_S(1,3),angu_pow_P(3,3)
      integer *4 angu_pow_D(6,3),angu_pow_F(10,3)
      integer *4 angu_pow_G(15,3)

      integer *4 ang_pow(ISPC,3)
      integer *4 n_pg_atom(ISPC)
      integer *4 k_m(3),k_n(3)
      integer *4 basis_n(ISPC)
    
      
      real *8 re(100,3),expon,coe
      real *8 ord_expon(ISPC),ord_coe(ISPC)
      real *8 R_m(4),R_n(4)
      real *8 R_ms(ISPC,3)

      DATA BLNK /'E','S','P','D','F','G'/
      DATA shell_p /'PX','PY','PZ'/
      DATA shell_d /'XX','YY','ZZ','XY','XZ','YZ'/
      DATA shell_f /'XXX','YYY','ZZZ','XXY','XXZ','YYX','YYZ',
     *     'ZZX','ZZY','XYZ'/
      DATA shell_G /'XXXX','YYYY','ZZZZ','XXXY','XXXZ','YYYX','YYYZ',
     *     'ZZZX','ZZZY','XXYY', 'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/

      
      DATA angu_pow_S /0,0,0/
      DATA angu_pow_P /1,0,0,
     *                 0,1,0,
     *                 0,0,1/
      DATA angu_pow_D /2,0,0,1,1,0,
     *                 0,2,0,1,0,1,
     *                 0,0,2,0,1,1/
      DATA angu_pow_F /3,0,0,2,2,1,0,1,0,1,
     *                 0,3,0,1,0,2,2,0,1,1,
     *                 0,0,3,0,1,0,1,2,2,1/
      DATA angu_pow_G /4,0,0,3,3,1,0,1,0,2,2,0,2,1,1,
     *                 0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,
     *                 0,0,4,0,1,0,1,3,3,0,2,2,1,1,2 /


      OPEN(unit=14,file='xyz.txt',STATUS='OLD',FORM='FORMATTED')
      OPEN(unit=15,file='capint.txt',STATUS='OLD',FORM='FORMATTED')
      
      DO i=1,ISPC
         Basis_N(i)=0
      END DO
      

      READ(14,*) n_atoms

      NG=1
      NSH=0
      iiii=1

      DO  k=1,n_atoms
         READ(14,*)IA,UW,(re(k,i),i=1,3)
         write(*,*)IA,UW,(re(k,i),i=1,3)
!         READ(15,*)(re(k,i),i=1,3)
         READ(15,*)IA
         write(*,*)IA
         Do ISNEW=1,200  
            READ(15,*)shell_type,n_con_guass
            write(*,*)shell_type,n_con_guass
            if(shell_type.eq.BLNK(1)) go to 201
            if(shell_type.eq.blnk(2)) go to 200
            if(shell_type.eq.blnk(3)) go to 300
            if(shell_type.eq.blnk(4)) go to 400
            if(shell_type.eq.blnk(5)) go to 500
            if(shell_type.eq.blnk(6)) go to 600
            
C!!   S TYPE ORBITALS 
            
 200        continue
            
       DO n_c_g=1,n_con_guass
          NSH=NSH+1
          READ(15,*)nj,expon,coe
          ord_expon(NSH)=expon
          ord_coe(NSH)=coe
          shell_arr(NSH)=shell_type
          do ia_p=1,3
             ang_pow(NSH,ia_p)=angu_pow_S(1,ia_p)
          end do
          basis_n(NG)=n_con_guass
       END DO
       NG=NG+1
       
       go to 1500
       
C!!   P TYPE ORBITALS
       
 300   continue
       
       DO n_c_g=1,n_con_guass
          NSH=NSH+1
          READ(15,*)nj,expon,coe
          do i=1,3
             shell_arr(NSH+(i-1)*n_con_guass)=shell_p(i)
             ord_expon(NSH+(i-1)*n_con_guass)=expon
             ord_coe(NSH+(i-1)*n_con_guass)=coe
             do ia_p=1,3
                ang_pow(NSH+(i-1)*n_con_guass,ia_p)=angu_pow_p(i,ia_p)
             end do
          end do                              
       END DO
       
       do  iji=1,3
          basis_n(NG)=basis_n(NG)+n_con_guass
          NG=NG+1
       end do
       
       NSH=NSH+2*n_con_guass
       go to 1500
       
C!!   D TYPE ORBITALS
       
 400   continue
       
       DO   n_c_g=1,n_con_guass
          NSH=NSH+1
          READ(15,*)nj,expon,coe
          do i=1,6
             shell_arr(NSH+(i-1)*n_con_guass)=shell_d(i)
             ord_expon(NSH+(i-1)*n_con_guass)=expon
             ord_coe(NSH+(i-1)*n_con_guass)=coe
             do ia_p=1,3
                ang_pow(NSH+(i-1)*n_con_guass,ia_p)=angu_pow_d(i,ia_p)
             end do
          end do
       END DO
       
       do  iji=1,6
          basis_n(NG)=basis_n(NG)+n_con_guass
          NG=NG+1
       end do
       NSH=NSH+5*n_con_guass
       go to 1500
       
C!!   F TYPE ORBITALS
 500   continue 
       
       DO   n_c_g=1,n_con_guass
          NSH=NSH+1
          READ(15,*)nj,expon,coe
          do i=1,10
             shell_arr(NSH+(i-1)*n_con_guass)=shell_f(i)
             ord_expon(NSH+(i-1)*n_con_guass)=expon
             ord_coe(NSH+(i-1)*n_con_guass)=coe
             do ia_p=1,3
                ang_pow(NSH+(i-1)*n_con_guass,ia_p)=angu_pow_f(i,ia_p)
             end do
          end do
       END DO
       
       do  iji=1,10
          basis_n(NG)=basis_n(NG)+n_con_guass
          NG=NG+1
       end do
       NSH=NSH+9*n_con_guass
       go to 1500


 600   continue
      DO   n_c_g=1,n_con_guass
          NSH=NSH+1
          READ(15,*)nj,expon,coe
          do i=1,15
             shell_arr(NSH+(i-1)*n_con_guass)=shell_f(i)
             ord_expon(NSH+(i-1)*n_con_guass)=expon
             ord_coe(NSH+(i-1)*n_con_guass)=coe
             do ia_p=1,3
                ang_pow(NSH+(i-1)*n_con_guass,ia_p)=angu_pow_g(i,ia_p)
             end do
          end do
       END DO
       
       do  iji=1,15
          basis_n(NG)=basis_n(NG)+n_con_guass
          NG=NG+1
       end do
       NSH=NSH+14*n_con_guass
       go to 1500



 1500   continue
       
      END DO 

 201       continue
 

C!!!MIRO-L1        
          n_pg_atom(k)=NSH

          If(k.eq.1)then
             INPG=0
          else 
             INPG= n_pg_atom(k-1)
          endif

          do ii= INPG+1,NSH

             do isw=1,3
                R_ms(ii,isw)=re(k,isw)
             end do
             iiii=iiii+1
          end do
       END DO
C!!!!MIRO-U1

       n_b=NG-1
       
C     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C       basis set  details for debugging
C       go to 602
       
       write(16,*)'number  of primitive guassians=',NSH
       write(16,*)'number of basis',n_b 

       KS=0
       DO i=1,N_B
          write(16,*)i,basis_n(i)
          do j=1,basis_n(i)
             KS=KS+1
             write(16,*)(R_MS(KS,l),l=1,3)
             write(16,*)KS,SHELL_ARR(KS),ORD_EXPON(KS),ORD_COE(KS)
           
          end do
       END DO

 602    continue

C!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~       

C!NORMALIZATION PART for Primitive Guassians
C!Only the primitives are normalized 
C!This normalization is compatible with the  
C!NORMF=1 NORMP=0 of gamess US. In this case the the SUBROUTINE
C!CONGUANORM  IS NOT NEEDED


        do i=1,NSH
           do isw=1,3
              R_m(isw)=R_ms(i,isw)
              k_m(isw)=(ang_pow(i,isw))
           end do
           
           R_m(4)=ord_expon(i)
           
           CALL PEM(R_m,R_m,K_m,K_m,PEINT)
           ord_coe(i)=(ord_coe(i)/(sqrt(PEINT)))
C           write (36,*)i,shell_arr(i),ord_expon(i),ord_coe(i)
C           write(16,*) (ang_pow(i,j),j=1,3)
C           write(16,*) (R_ms(i,j),j=1,3)
        end do
        RETURN
        END

C*******************************************************************************
      SUBROUTINE CONGUANORM(N_B,ISPC,ord_expon,ord_coe,R_ms,basis_n,
     *     ang_pow) 
      
      IMPLICIT REAL*8 (A-H,O-Z)

      integer *4 n_b
      integer *4 ang_pow(ISPC,3)
      integer *4 k_m(3),k_n(3)
      integer *4 basis_n(ISPC)
      integer *4 basis_ord(ISPC)
     
      real *8 ord_expon(ISPC),ord_coe(ISPC)
      real *8 R_m(4),R_n(4)
      real *8 R_ms(ISPC,3)
      real *8 PW
      real *8 PE_AO,PEINT
      
      EXTERNAL IPOINT

C! THis is a subroutine for Normalizing the Contracted guassians.

      do i=1,n_b
         if(i.eq.1)then
            iasp=0
         else
            iasp=basis_ord(i-1)
         endif
         basis_ord(i)=iasp+basis_n(i)
      end do
      
CMI .... run over contracted basis functions ...

      DO ii=1,n_b
         jj=ii
        
         IPU=IPOINT(ii,jj)
         PE_AO=0.0d0
         if (ii.eq.1)then
            iangii=0
         else
            iangii=basis_ord(ii-1)
         endif
         if (jj.eq.1)then
            iangjj=0
         else
            iangjj=basis_ord(jj-1)
         endif
         
         Do i=iangii+1,basis_ord(ii)
            do j=iangjj+1,basis_ord(jj)
               
               do isw=1,3
                  R_m(isw)=R_ms(i,isw)
                  R_n(isw)=R_ms(j,isw)
                  k_m(isw)=(ang_pow(i,isw))
                  k_n(isw)=(ang_pow(j,isw))
               end do

               R_m(4)=ord_expon(i)
               R_n(4)=ord_expon(j)
               
               COMUL=ord_coe(i)*ord_coe(j)
               CALL PEM(R_m,R_n,K_m,K_n,PEINT)
               PE_AO=PE_AO+COMUL*PEINT

            end do
         End Do

         Do i=iangii+1,basis_ord(ii)
            unnormcoe=ord_coe(i)
            ord_coe(i)=unnormcoe/sqrt(PE_AO)
         EndDo
         
      END DO
      
      RETURN
      END
     
C ===============================================================
      SUBROUTINE SMATRIX(N_B,ISPC,ord_expon,ord_coe,R_ms,basis_n,
     *     ang_pow,PE_AO) 

      IMPLICIT REAL*8 (A-H,O-Z)

      integer *4 n_b
      integer *4 ang_pow(ISPC,3)
      integer *4 k_m(3),k_n(3)
      integer *4 basis_n(ISPC)
      integer *4 basis_ord(ISPC)
     
      real *8 ord_expon(ISPC),ord_coe(ISPC)
      real *8 R_m(4),R_n(4)
      real *8 R_ms(ISPC,3)
      real *8 PW

      real *8 PEINT,PE_AO(N_B*(N_B+1)/2)
      
      EXTERNAL IPOINT
      
C!THis subroutine generates the overlap matrix in AO basis.
C!THis is to compare the nomalization conditions with the 
C!commercial codes or the codes written by a 3rd party.


      do i=1,n_b
         do j=1,n_b
            IPU=IPOINT(i,j)
            PE_AO(IPU)=0.d0
         end do
      end do

      do i=1,n_b
         if(i.eq.1)then
            iasp=0
         else
            iasp=basis_ord(i-1)
         endif
         basis_ord(i)=iasp+basis_n(i)
      end do
      
      DO ii=1,n_b

         DO jj=1,ii

            IPU=IPOINT(ii,jj)
            
            if (ii.eq.1)then
               iangii=0
            else
               iangii=basis_ord(ii-1)
            endif
            if (jj.eq.1)then
               iangjj=0
            else
               iangjj=basis_ord(jj-1)
            endif
            
            Do i=iangii+1,basis_ord(ii)
               do j=iangjj+1,basis_ord(jj)
                  
                  do isw=1,3
                     R_m(isw)=R_ms(i,isw)
                     R_n(isw)=R_ms(j,isw)
                     k_m(isw)=(ang_pow(i,isw))
                     k_n(isw)=(ang_pow(j,isw))
                  end do
                  R_m(4)=ord_expon(i)
                  R_n(4)=ord_expon(j)
                  
                  COMUL=ord_coe(i)*ord_coe(j)
C                  write(*,*)ord_coe(i),ord_coe(j)
                  CALL PEM(R_m,R_n,K_m,K_n,PEINT)
                  PE_AO(IPU)= PE_AO(IPU)+COMUL*PEINT
               end do
            End Do
         END DO
      END DO
      
      open(unit=14,file='SMATRIX',STATUS='UNKNOWN',FORM='FORMATTED')

      do ii=1,n_b
         do jj=1,ii
            IPU=IPOINT(ii,jj)
            if (jj.gt.ii)go to 301
            write(14,*)ii,jj,PE_AO(IPU)
         end do
 301     continue
      end do

      RETURN
      END
     

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION IPOINT(ii,jj)
      IMPLICIT REAL*8 (A-H,O-Z)
      if(ii.ge.jj)then
         IPOINT=(ii*(ii-1)/2)+jj
      else
         IPOINT=(JJ*(JJ-1)/2)+ii
      endif
      END 

C**********************************************

      subroutine SMAT(A,B,L1,L2,OVINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer *4 L1, L2
      real *8  A(2),B(2)
      real *8 GAMMA,P,PA,PB, AB, ABSQ
      real *8 OVINT,RMULT,PREMULT, PI
      INTEGER *4 IFI
      real *8 GAMMA1, GAMMA2
      real *8 BINC
      data PI/3.141592653589793/
      EXTERNAL  ISUM
   
      GAMMA=(A(2))+(B(2))
      P=((A(1))*(A(2))+(B(1))*(B(2)))/GAMMA
      PA=P-A(1)
      PB=P-B(1)
      AB=(A(1))-(B(1)) 
      OVINT=0.0d0

      DO i=0,(L1+l2)/2
         IFI=ISUM(2*i)
         GAMMA1=(2.0D0*GAMMA)**(i)
         GAMMA2=(PI/GAMMA)**(0.50d0)
         RMULT=IFI*GAMMA2/GAMMA1
         CALL BINCOEF(2*i,L1,L2,PA,PB,BINC)
         OVINT=OVINT+(RMULT*BINC) 
      END DO
      ABSQ=AB*AB
      PREMULT=exp((-(A(2))*(B(2))*ABSQ)/gamma)
      OVINT=OVINT*PREMULT

      RETURN
      END


C*******************************************************

      SUBROUTINE PEM(A,B,KA,KB,PEINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      
      real *8 A(4),B(4)
      integer *4 KA(3),KB(3)
      real *8 PEINT,PREMULT
      real *8 GAMMA
      real *8 AB(3),P(3),PA(3),PB(3)
      REAL *8 OVINT(3)
      
      GAMMA=(A(4))+(B(4))
      DO ig=1,3
         AB(ig)=(A(ig))-(B(ig))
         P(ig)=((A(4))*(A(ig)) +(B(4))*(B(IG)))/GAMMA
         PA(IG)=((P(IG))-(A(IG)))
         PB(IG)=((P(IG))-(B(IG)))
      END DO
      
      ABSQ=(AB(1))*(AB(1))+(AB(2))*(AB(2))+(AB(3))*(AB(3))
      PREMULT=DEXP(((-A(4))*(B(4))*ABSQ)/GAMMA)
      
      DO I=1,3
         KAI=KA(I)
         KBI=KB(I)
         PAI=PA(I)
         PBI=PB(I)
         CALL OVMAT(GAMMA,KAI,KBI,PAI,PBI,OVINTI)
         OVINT(I)=OVINTI
      END DO
      
      PEINT=PREMULT*(OVINT(1))*(OVINT(2))*(OVINT(3))

      RETURN
      END


      SUBROUTINE OVMAT(GAMMA,L1,L2,AP,BP,OVINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 LI,L2
      REAL *8 GAMMA,KA,KB,AP
      REAL *8 PI,GAMMA1,GAMMA2,PREMULT
      REAL *8 BINC
      INTEGER *4 IFI
      EXTERNAL ISUM
      data PI/3.141592653589793/

      OVINT=0.0d0
      DO i=0,(L1+L2)/2
         IFI=ISUM(2*i)
         GAMMA1=(2.0D0*GAMMA)**(i)
         GAMMA2=(PI/GAMMA)**(0.50d0)
         PREMULT=IFI*GAMMA2/GAMMA1
         CALL BINCOEF(2*i,L1,L2,AP,BP,BINC)
         OVINT=(OVINT)+(PREMULT*BINC)
      END DO
      RETURN
      END

      SUBROUTINE  BINCOEF(K,L1,L2,AP,BP,BINC)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 K,L1,L2
      REAL *8 AP,BP,BINC
      REAL *8 COEL1,COEL2,PRM1,PRM2,PRM3
      REAL *8 PREMULT
      EXTERNAL BINEXP
      BINC=0.0D0

      DO I=0,L1
         DO J=0,L2
            IF((I+J).EQ.K)THEN
               COEL1=BINEXP(L1,I)
               COEL2=BINEXP(L2,J)
               PRM1=COEL1*COEL2
               IL=L1-i
               JL=L2-J
               PRM2=(AP)**(IL)
               PRM3=(BP)**(JL)
               PREMULT=PRM1*PRM2*PRM3
               BINC=BINC+PREMULT
            ENDIF
         END DO
      END DO
      RETURN
      END


      FUNCTION BINEXP(L,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 BINEXP
      INTEGER *4 L,I,LP1,LP2,L3
      INTEGER *4 NFACT
      EXTERNAL NFACT
      LP1=NFACT(L)
      LP2=NFACT(I)
      LP3=NFACT(L-I)
      BINEXP=(LP1)/(LP2*LP3)
      END

      FUNCTION NFACT(I)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 I,K,KL,NFACT
      KL=1
      DO K=1,I
         KL=KL*K
      END DO
      NFACT=KL
      END

      FUNCTION ISUM(I)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 I,K,KL,M,ISUM
      KL=1
      Do k=1,i
         m=(k-((k/2)*2))
         if(m.eq.0)go to 110
         KL=kL*k
 110     continue
      End Do
      ISUM=KL
      END
      
      SUBROUTINE MOINTEGRAL(N_B,AOMO,AOINT,MOINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      real *8 AOMO(N_B,N_B)
      real *8 AOINT((N_B*(N_B+1))/2)
      real *8 MOINT((N_B*(N_B+1)/2))
      EXTERNAL IPOINT
      
      Do i=1,n_b
         do j=1,i
            IPU=IPOINT(i,j)
            MOINT(IPU)=0.d0
         end do
      End do
      
      DO ii=1,n_b
         Do jj=1,ii
            IPU=IPOINT(ii,jj)
            do i=1,n_b
               do j=1,n_b
                  IPS=IPOINT(i,j)
                  if(j.gt.i)IPS=IPOINT(j,i)
                  COMULT=AOMO(jj,j)*AOMO(ii,i)
                  MOINT(IPU)= MOINT(IPU)+COMULT*AOINT(IPS)
               end do
            end do
         End  Do
      END DO
      
      RETURN
      END

!N - number of AO, N1 - number of MO
      SUBROUTINE CITR(N,N1,C,H,CHAUX)
      IMPLICIT REAL*8 (A-H,O-Z)
      real *8 B(N,N), B1(N1,N1)
      real *8 C(N1,N),D(N1,N)
      real *8 H((N*(N+1))/2)

      character *10  CHAUX
      EXTERNAL IPOINT
      CHARACTER *1 TRANSA, TRANSB
      alpha=1.0d0
      beta=0.0d0

      do i=1,N
      do j=1,N
         IPU=IPOINT(i,j)
         B(i,j)=H(IPU)
      end do
      end do

      CALL DGEMM('N','N',N1,N,N,alpha,C,N1,B,N,beta,D,N1)
      CALL DGEMM('N','T',N1,N1,N,alpha,D,N1,C,N1,beta,B1,N1)
      open(unit=18,file=CHAUX,STATUS='UNKNOWN',FORM='FORMATTED')

      write(18,*)n1
      do ii=1,n1
         do jj=1,ii
            IPU=IPOINT(ii,jj)
            write(18,211)B1(ii,JJ)
            write(86,*)ii,jj,B1(ii,JJ)
 211        format( f25.14)
         end do
      end do  

      RETURN
      END


      SUBROUTINE AOINTEGRAL(N_B,ISPC,ord_expon,ord_coe,R_ms,basis_n,
     *     ang_pow,PE_AO)      

      
      
      IMPLICIT REAL*8 (A-H,O-Z)

      integer *4 n_b
      integer *4 ang_pow(ISPC,3)
      integer *4 k_m(3),k_n(3)
      integer *4 basis_n(ISPC)
      integer *4 basis_ord(ISPC)
     
      real *8 ord_expon(ISPC),ord_coe(ISPC)
      real *8 R_m(4),R_n(4)
      real *8 R_ms(ISPC,3)
          
      
      real *8 PEINT,PE_AO((N_B*(N_B+1))/2)
      real *8 core(3),cen(3)
      common /boxsize/ core
      common /symcent/cen
      EXTERNAL IPOINT
     
     
      do i=1,n_b
         do j=1,i
            IPU=IPOINT(i,j)
            PE_AO(IPU)=0.d0
         end do
      end do

      do i=1,n_b
         if(i.eq.1)then
            iasp=0
         else
            iasp=basis_ord(i-1)
         endif
         basis_ord(i)=iasp+basis_n(i)
      end do

      DO ii=1,n_b
C     write(*,*)ii
         DO jj=1,ii
            IPU=IPOINT(ii,jj)
            
            if (ii.eq.1)then
               iangii=0
            else
               iangii=basis_ord(ii-1)
            endif
            if (jj.eq.1)then
               iangjj=0
            else
               iangjj=basis_ord(jj-1)
            endif
            
            Do i=iangii+1,basis_ord(ii)
               do j=iangjj+1,basis_ord(jj)
                  
                  do isw=1,3
                     R_m(isw)=R_ms(i,isw)
                     R_n(isw)=R_ms(j,isw)
                     k_m(isw)=(ang_pow(i,isw))
                     k_n(isw)=(ang_pow(j,isw))
                  end do
                  R_m(4)=ord_expon(i)
                  R_n(4)=ord_expon(j)
                  
                  COMUL=ord_coe(i)*ord_coe(j)
C                  write(*,*)ord_coe(i),ord_coe(j)
                  CALL CAPM(R_m,R_n,K_m,K_n,PEINT)
                  PE_AO(IPU)= PE_AO(IPU)+COMUL*PEINT
C                  write(*,*)PEINT,1/COMUL
               end do
            End Do
         END DO
      END DO
      
    
C      write(*,*)'integrals in AO basis'

      RETURN
      END

 
