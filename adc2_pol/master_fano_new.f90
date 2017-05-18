subroutine master_fano()

  use constants
  use parameters
  use select_fano
!!$  use get_matrix
  use fspace
  use partgammas
  use misc
  use davmod
  
  implicit none

  integer, dimension(:,:), allocatable :: kpq_init,kpq_fin
  integer, dimension(:), allocatable :: inddx
  integer :: i,ndimi,ndimf,nstate,nstatef,nout,k,ii,nosstates,real_lancstates,ccnt
  integer*4 :: noffdelem

  real(d) :: sum,pgamma,cntr,matt

  real(d), dimension(:), allocatable :: ener,isen,fsen,mtm,tmvec,nfsen,fsens,ffsen,gammavec,mgammavec,fgammavec,tempvecc
  real(d), dimension(:,:), allocatable :: arr,isfamani,fsfamani

  external M01CBF

  if((fmethod .lt. 1) .or. (fmethod .gt. 3)) then
     write(6,*) "You have not chosen Fano space. Set the field *fanom to 1, 2, or 3."
     stop
  end if

  select case(fmethod)

!!$*****************************************************************************************************
!!$***************************************FANO TDA******************************************************
!!$*****************************************************************************************************

  case(1)
     
     write(6,*) "Constructing Fano space in TDA approximation"
     
     allocate(kpq_init(7,0:nBas**2*4*nOcc**2),kpq_fin(7,0:nBas**2*4*nOcc**2)) 
     
     call select_singles(kpq_init(:,:),-1)
     call select_singles(kpq_fin(:,:),1)
     
     ndimi=kpq_init(1,0)
     ndimf=kpq_fin(1,0)
     
     write(6,*) 'Initial space dim',ndimi
     write(6,*) 'Final space dim',ndimf

     if (chrun .eq. 'dire') then
        
        write(6,*) 'Direct diagonalization of the initial and final subspaces'
        
        allocate(isfamani(ndimi,numinista),isen(numinista),arr(ndimi,ndimi),ener(ndimi))
        call get_fstates_direct(ndimi,kpq_init(:,:),arr(:,:),ener(:),nstate,'tda','i')
        isfamani(:,:)=arr(:,1:numinista)
        isen(:)=ener(1:numinista)
        deallocate(arr,ener)
        allocate(arr(ndimf,ndimf),ener(ndimf))
        call get_fstates_direct(ndimf,kpq_fin(:,:),arr(:,:),ener(:),nstate,'tda','f')
        allocate(fsfamani(ndimf,nstate),fsen(nstate))
        fsfamani(:,:)=arr(:,1:nstate)
        fsen(:)=ener(1:nstate)
        deallocate(arr,ener)

        allocate(gammavec(nstate))
        call getf_HIJ_tda(ndimi,ndimf,numinista,nstate,kpq_init(:,:),kpq_fin(:,:),&
             isfamani(:,:),isen(:),fsfamani(:,:),gammavec(:))
!        call get_gamma_1(nstate,numinista,isen(:),fsen(:),gammavec(:))
        deallocate(gammavec,isfamani,isen,fsfamani,fsen)
        
        deallocate(kpq_init,kpq_fin)
        
     elseif (chrun .eq. 'save') then

        write(6,*) 'Only full diagonalisation in the TDA case'
        stop
        
     end if
        
     
 

!!$***************************************************************************************************************
!!$*****************************************FANO ADC2*************************************************************
!!$***************************************************************************************************************
  case(2)
     write(6,*) "Constructing Fano space in ADC2 approximation"

     if (info .eq. 1) then
        allocate(kpq_init(7,0:nBas**2*4*nOcc**2))
        call  select_atom_is(kpq_init(:,:))
        call  select_atom_d(kpq_init(:,:),-1)
        ndimi=kpq_init(1,0)+kpq_init(2,0)+kpq_init(3,0)+kpq_init(4,0)+2*kpq_init(5,0)
        write(6,*) 'Initial space dim',ndimi
        deallocate(kpq_init)
        allocate(kpq_fin(7,0:nBas**2*4*nOcc**2))
        call  select_atom_fs(kpq_fin(:,:))
        call  select_atom_d(kpq_fin(:,:),1)
        ndimf=kpq_fin(1,0)+kpq_fin(2,0)+kpq_fin(3,0)+kpq_fin(4,0)+2*kpq_fin(5,0)
        write(6,*) 'Final space dim',ndimf
        deallocate(kpq_fin)
        stop
     end if


     if (chrun.eq.'direct') then
        
        write(6,*) 'ADC2 Direct diagonalisation is not available on ICD cluster'
        stop

!============== ADC2 INI ============================
     elseif (chrun .eq. 'save') then

        allocate(kpq_init(7,0:nBas**2*4*nOcc**2))
        call  select_atom_is(kpq_init(:,:))
        call  select_atom_d(kpq_init(:,:),-1)
        ndimi=kpq_init(1,0)+kpq_init(2,0)+kpq_init(3,0)+kpq_init(4,0)+2*kpq_init(5,0)
        write(6,*) 'Initial space dim',ndimi 

        if (idiag .eq. 1) then
          write(6,*) 'ADC2 IS Full diagonalisation is not available on ICD cluster'
          stop  
        elseif (idiag .eq. 2) then

           call write_fspace_adc2_1(ndimi,kpq_init(:,:),noffdelem,'i')
           allocate(isfamani(ndimi,davstates),isen(davstates))
           write(*,*) 'davstates',davstates
           
           dmain=kpq_init(1,0)
           mtxidd='init'
           call master_dav(ndimi,noffdelem,'i')
           call get_bound(ndimi,davstates,kpq_init(:,:),isfamani(:,:),isen(:),nstate)

           if (tranflag .eq. 'y') then
             allocate(mtm(ndimi),tmvec(nstate))
             write(6,*) 'Calculating transition moments in ',tranmom,' direction. for the initial manifold'
             call get_modifiedtm_adc2(ndimi,kpq_init(:,:),mtm(:))
              do i=1,nstate
                 tmvec(i)=dot_product(isfamani(:,i),mtm(:))
              end do
              write(6,*) 'The following initial states were selected'
              call table2(ndimi,nstate,isen(1:nstate),isfamani(:,1:nstate),tmvec(:),tmvec(:))
              deallocate(mtm,tmvec)
           end if

         write(6,*) 'The requested Fano state lies at energy of',isen(ninista),' a.u.'

        else 
           write(6,*) 'Set idiag to 1 or 2'
           stop
        end if

!============ ADC2 FINAL ============================= 
        allocate(kpq_fin(7,0:nBas**2*4*nOcc**2))
        call  select_atom_fs(kpq_fin(:,:))
        call  select_atom_d(kpq_fin(:,:),1)
        ndimf=kpq_fin(1,0)+kpq_fin(2,0)+kpq_fin(3,0)+kpq_fin(4,0)+2*kpq_fin(5,0)
        write(6,*) 'Final space dim',ndimf
        
        if (fdiag .eq. 1) then

!<<<<<<<<<<<<<<<<< FS Full Diag >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

         allocate(fsfamani(ndimf,ndimf),fsen(ndimf))
         call get_fstates_direct(ndimf,kpq_fin(:,:),fsfamani(:,:),fsen(:),nstatef,'ad2','f')
    
         allocate(mgammavec(ndimf))
         call getf_HIJ_adc2(ndimi,ndimf,ninista,nstate,kpq_init,kpq_fin,isen,isfamani,mgammavec)
         deallocate(isfamani,isen)
       
         allocate(gammavec(ndimf),tempvecc(ndimf),nfsen(ndimf))
         ccnt=0

         do i=1,ndimf
              do ii=1,kpq_fin(1,0)
                cntr=cntr+fsfamani(ii,i)**2
              end do
             if (fsen(i).lt.4) then
             if (cntr.gt.0.05) then
             ccnt=ccnt+1
             tempvecc(:)=fsfamani(:,i)
             nfsen(ccnt)=fsen(i)
             pgamma=dsp(ndimf,mgammavec(:),tempvecc(:))
             gammavec(ccnt)=2._d*pi*pgamma**2
             end if
             end if
             cntr=0.0
         end do

         deallocate(fsfamani,mgammavec,tempvecc,fsen)

         call get_gamma_1(ccnt,nfsen,gammavec)

         deallocate(nfsen,gammavec)

         write(*,*) 'ADC2 FANO: Save d(I):2 - d(F):1 ended successfully'

!<<<<<<<<<<<<<<<< FS Lanc Diag >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      elseif (fdiag .eq. 2) then

           call write_fspace_adc2_1(ndimf,kpq_fin(:,:),noffdelem,'f')

           write(*,*)'-----ADC2  Info-----'
           write(*,*) 'ndimf',ndimf
           write(*,*) 'ndimi',ndimi
           write(*,*) 'nstate',nstate
           write(*,*) 'ninista',ninista
           write(*,*) 'noffdiag',noffdelem
           write(*,*)'--------------------'

           allocate(mgammavec(ndimf))
           call getf_HIJ_adc2(ndimi,ndimf,ninista,nstate,kpq_init,kpq_fin,isen,isfamani,mgammavec)
           nout=kpq_fin(1,0)
           deallocate(kpq_init,kpq_fin,isfamani)
           write(*,*) 'ADC2 GammaVec built successfully'

           if(lmain.gt.nout) then
           write(*,*) 'lmain must not be greater than nout!'
           stop
           end if
           if (lmain.gt.1000) then
           write(*,*) 'lmain is too large!' !can be changed in parameter.f90
           stop
           else
           write(*,*) 'lmain size is valid'
           end if
           call fill_stvc(nout,mgammavec(1:nout))
           mtxidl='fina'
           call master_lancdiag(ndimf,noffdelem,'f')
           !stop   ! ============TEST!!!!!!!!!!!!!!!!!!!!!!!         

          
           allocate(fsen(lancstates))
           call readen(ndimf,2,lancstates,fsen,real_lancstates)
           if(real_lancstates.lt.lancstates) then
           write(*,*) 'Attention: Not all lancstates available!'
           end if
           
           allocate(gammavec(real_lancstates),nfsen(real_lancstates))
           call read_gammavec(ndimf,lancstates,real_lancstates,nout,mgammavec(:),fsen(:),nosstates,nfsen(:),gammavec(:))
           deallocate(mgammavec,fsen)

!Anpassung von nfsen und gammavec
              allocate(fgammavec(nosstates),ffsen(nosstates))
              do ii=1,nosstates
               ffsen(ii)=nfsen(ii)
               fgammavec(ii)=gammavec(ii)
              end do
              deallocate(nfsen,gammavec)
!Anpassung abgeschlossen

          write(*,*) 'Partgammas at the energy',isen(ninista)
          call get_gamma_1(nosstates,ffsen(:),fgammavec(:))
          deallocate(fgammavec,isen,ffsen)

        end if        
        
     end if
     
!!$*****************************************************************************************************
!!$****************************************FANO ADC2e***************************************************
!!$*****************************************************************************************************

  case(3)
     
     write(6,*) "Constructing Fano space in ADC2-ext.approximation"

     if (info .eq. 1) then
        allocate(kpq_init(7,0:nBas**2*4*nOcc**2))
        call select_atom_is(kpq_init(:,:))
        call select_atom_d(kpq_init(:,:),-1)
        ndimi=kpq_init(1,0)+kpq_init(2,0)+kpq_init(3,0)+kpq_init(4,0)+2*kpq_init(5,0)
        write(6,*) 'Initial space dim',ndimi
        deallocate(kpq_init)
        allocate(kpq_fin(7,0:nBas**2*4*nOcc**2))

        call  select_atom_fs(kpq_fin(:,:))
        call  select_atom_d(kpq_fin(:,:),1)
        ndimf=kpq_fin(1,0)+kpq_fin(2,0)+kpq_fin(3,0)+kpq_fin(4,0)+2*kpq_fin(5,0)
        write(6,*) 'Final space dim',ndimf
        deallocate(kpq_fin)

!        write(*,*) '---ANA----'
!        write(*,*) '1:12(1-96)'
!        write(*,*) '.1:12(3-121)'
!        write(*,*) 'D1.1',vpqrs(1,96,3,121)
!        write(*,*) 'E1.1',vpqrs(1,3,96,121)
!        write(*,*) '.2:12(3-122)'
!        write(*,*) 'D1.1',vpqrs(1,96,3,122)
!        write(*,*) 'E1.1',vpqrs(1,3,96,122)
!        write(*,*) '.3:12(3-102)'
!        write(*,*) 'D1.1',vpqrs(1,96,3,102)
!        write(*,*) 'E1.1',vpqrs(1,3,96,102)
!         matt=0.0
!         call adc2ext_0_0(21,96,1,121,3,matt)
!         write(*,*) 'Mat 1.1',matt
!         matt=0.0
!         call adc2ext_0_0(21,96,1,122,3,matt)
!         write(*,*) 'Mat 1.2',matt
        stop
     end if
     
     if (chrun .eq. 'dire') then
        
      write(6,*) 'ADC2e Direct diagonalisation not available!'
      stop      
  
     elseif (chrun .eq. 'save') then

!==============INI ADC2e======================================
        allocate(kpq_init(7,0:nBas**2*4*nOcc**2))
        call select_atom_is(kpq_init(:,:))
        call select_atom_d(kpq_init(:,:),-1)
        ndimi=kpq_init(1,0)+kpq_init(2,0)+kpq_init(3,0)+kpq_init(4,0)+2*kpq_init(5,0)
        write(6,*) 'Initial space dim',ndimi 

        if (idiag .eq. 1) then
         
          write(6,*) 'ADC2e Full IS diagonalisation not available!'
          stop

        elseif (idiag .eq. 2) then
           call write_fspace_adc2e_1(ndimi,kpq_init(:,:),noffdelem,'i')
           allocate(isfamani(ndimi,davstates),isen(davstates))
            
           dmain=kpq_init(1,0)
           mtxidd='init'
           call master_dav(ndimi,noffdelem,'i')
           call get_bound(ndimi,davstates,kpq_init(:,:),isfamani(:,:),isen(:),nstate)

           if (tranflag .eq. 'y') then
              allocate(mtm(ndimi),tmvec(nstate))
              write(6,*) 'Calculating transition moments in ',tranmom,' direction. for the initial manifold'
              call get_modifiedtm_adc2(ndimi,kpq_init(:,:),mtm(:))
              do i=1,nstate
                 tmvec(i)=dot_product(isfamani(:,i),mtm(:))
              end do
              write(6,*) 'The following initial states were selected'
              call table2(ndimi,nstate,isen(1:nstate),isfamani(:,1:nstate),tmvec(:),tmvec(:))
              deallocate(mtm,tmvec)
           end if

           write(6,*) 'The requested Fano state lies at energy of',isen(ninista),' a.u.'

        else 
           write(6,*) 'Set idiag to 1 or 2'
           stop
        end if

!--------------------------Final State (ADC2e, save)

        !Obtaining final states 
        allocate(kpq_fin(7,0:nBas**2*4*nOcc**2))
        call  select_atom_fs(kpq_fin(:,:))
        call  select_atom_d(kpq_fin(:,:),1)
        ndimf=kpq_fin(1,0)+kpq_fin(2,0)+kpq_fin(3,0)+kpq_fin(4,0)+2*kpq_fin(5,0)
        write(6,*) 'Final space dim',ndimf
        
        if (fdiag .eq. 1) then
 
!==================== ADC2e FS Full Diag ================

         allocate(fsfamani(ndimf,ndimf),fsen(ndimf))
         call get_fstates_direct(ndimf,kpq_fin(:,:),fsfamani(:,:),fsen(:),nstatef,'a2e','f')

         allocate(mgammavec(ndimf))
!         call getf_HIJ_adc2(ndimi,ndimf,ninista,nstate,kpq_init,kpq_fin,isen,isfamani,mgammavec)
         call get_mgamma_adc2e(ndimi,ndimf,nstate,kpq_init,kpq_fin,isfamani,isen,mgammavec)
         deallocate(isfamani,isen)

         allocate(gammavec(ndimf),tempvecc(ndimf),nfsen(ndimf))
         ccnt=0
         do i=1,ndimf
              do ii=1,kpq_fin(1,0)
                cntr=cntr+fsfamani(ii,i)**2
              end do
             if (fsen(i).lt.4) then
             if (cntr.gt.0.05) then
             ccnt=ccnt+1
             tempvecc(:)=fsfamani(:,i)
             nfsen(ccnt)=fsen(i)
             pgamma=dsp(ndimf,mgammavec(:),tempvecc(:))
             gammavec(ccnt)=2._d*pi*pgamma**2
             end if
             end if
             cntr=0.0
         end do
         write(*,*) 'State Filter:',ccnt,'1h1p states found'
         deallocate(fsfamani,mgammavec,tempvecc,fsen)

         call get_gamma_1(ccnt,nfsen,gammavec)

         deallocate(nfsen,gammavec)

         write(*,*) 'ADC2ext FANO: fmeth: 3, Save d(I):2 - d(F):1 ended successfully'

!==================== ADC2e FS Lanc ====================
        elseif (fdiag .eq. 2) then
           call write_fspace_adc2e_1(ndimf,kpq_fin(:,:),noffdelem,'f')

           write(*,*)'-----ADC2e Info-----'
           write(*,*) 'ndimf',ndimf
           write(*,*) 'ndimi',ndimi
           write(*,*) 'nstate',nstate
           write(*,*) 'ninista',ninista
           write(*,*) 'noffdiag',noffdelem
           write(*,*)'--------------------'

           allocate(mgammavec(ndimf))
           call get_mgamma_adc2e(ndimi,ndimf,nstate,kpq_init,kpq_fin,isfamani,isen,mgammavec(:))
           nout=kpq_fin(1,0)
           deallocate(kpq_init,kpq_fin,isfamani)
           write(*,*) 'ADC2e GammaVec built successfully'

           if(lmain.gt.nout) then                            ! für die 1h1p variante/ pRICD
           write(*,*) 'lmain must not be greater than nout!'
           stop
           end if
           if (lmain.gt.1000) then
           write(*,*) 'lmain is too large!' !can be changed in parameter.f90
           stop
           else
           write(*,*) 'lmain size is valid'
           end if
           call fill_stvc(nout,mgammavec(1:nout))
           mtxidl='fina'
           call master_lancdiag(ndimf,noffdelem,'f')


           allocate(fsen(lancstates))
           call readen(ndimf,2,lancstates,fsen,real_lancstates)
           if(real_lancstates.lt.lancstates) then
           write(*,*) 'Attention: Not all lancstates available!'
           end if

           allocate(gammavec(real_lancstates),nfsen(real_lancstates))
           call read_gammavec(ndimf,lancstates,real_lancstates,nout,mgammavec(:),fsen(:),nosstates,nfsen(:),gammavec(:))
           deallocate(mgammavec,fsen)

!Anpassung von nfsen und gammavec
              allocate(fgammavec(nosstates),ffsen(nosstates))
              do ii=1,nosstates
               ffsen(ii)=nfsen(ii)
               fgammavec(ii)=gammavec(ii)
              end do
              deallocate(nfsen,gammavec)
!Anpassung abgeschlossen

            call get_sums(nosstates,ffsen(:),fgammavec(:))

            write(*,*) 'Partgammas at the energy',isen(ninista)
            call get_gamma_1(nosstates,ffsen(:),fgammavec(:))
            deallocate(fgammavec,isen,ffsen)

           end if
        end if
  end select
  
end subroutine master_fano
