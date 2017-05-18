      subroutine stieltjes_phi (epsiv,omega,num,e8,g8,gamma0,unit)
c Vitali Averbukh (2003). Comments to: vitali@tc.pci.uni-heidelberg.de
c This routine receives the sequence of the energy - oscillator strength pairs and 
c transforms it into a shorter sequence of the Stieltjes energy - oscillator strength
c pairs using the algorithm of Muller-Plathe & Dierksen [Electronic Structure of Atoms, 
c Molecules and Solids, Proceedings of the II Escola Brasileira de Estruture Eletronica 
c (World Scientific, Singapore, 1990), p.1-29]. The tutorial is available on the net: 
c                              http://citeseer.nj.nec.com/465113.html
c Sequences of more than three points and up to NP points are allowed. 
c The recursive relations for the coefficients of the orthogonal polynomials
c are implemented in quadruple precision. The maximal polynomial order is determined 
c by the magnitude of the polynomial overlap approacing the polynomial norm.
c If only low-order (MAXORD < 5) approximation is available, the program produces the 
c interpolated f(E) function. Otherwise, a sequence of f(E) functions for the 
c approximation orders from 5 to MAXORD so that the convergence can be seen.
c The tridiagonal coefficient matrix is diagonalized by the REAL*16 version of the TQL2 
c routine (source available).
c The cumulative F is differentiated numerically and the resulting set of points is 
c interpolated using monotonicity-preserving piecewise cubic Hermite interpolant 
c (NAG: E01BEF,E01BFF). This interpolation scheme ensures that the f(E) 
c is non-negative. The gamma at the ENERGY(1) is calculated by successively 
c higher order Stieltjes approximation. By default, the average of the three highest-order 
c results is taken as the final one. If convergence is detected at some lower order, the 
c converged result is preferred.
c The following parameters are transferred from the main program:
c OMEGA - the photon energy of interest
c CONV_MAX - the maximal number of convergence search loops
c OVERMAX -  the maximal permitted norm-to-overlap ratio of the 
c  adjacent orthogonal polynomials
c CONV_THRESH0 - initial convergence threshold
c CONV_FAC - convergence threshold factor making the convergence requirement
c more strict for lower Stieltjes orders
      implicit none

      integer*2 num
      integer np,maxord,ierr,imax,min,max,nint,unit
      integer k,n,i,j,ifail,converge,conv_max,iconv

c the maximal number of the energy points is NP
      parameter (np=3000,nint=200)

      real*16 e_point(np),g_point(np) 
      real*16 e_new(np),g_new(np)
      real*16 bcoef(0:np),acoef(np),qpol(0:np,np)
      real*16 diag(np),offdiag(np),abvec(np,np)
      real*16 e_min,e_max,bprod,asum,qnorm,qoverlap,overmax

      real*8 epsiv,omega,e8(np),g8(np),e0(np),gamma(np),der(np)
      real*8 energy(1),gammas(1),gamma0,gamma_h(np)
      real*8 energy_inter(nint+1),gammas_inter(nint+1)
      real*8 dif12,dif23,dif13,gmean,difmean
      real*8 conv_thresh,conv_thresh0,conv_fac

c      common /stieltjes_p/ overmax,conv_thresh0,
c     &     conv_fac,conv_max

c check the number of the input energy points 

      conv_max=10
      overmax=100.0
      conv_thresh0=0.05
      conv_fac=1.0
      
      if (num.gt.np) then
         print*, '***WARNING*** Stieltjes: too many energy points'
         print11, 'NUM (',num,') > ',np
         gamma0=0.d0
         return
      end if
 11   format(10x,a5,i3,a4,i3)
      if (num.le.3) then
         print*, '***ERROR*** Stieltjes: not enough energy points'
         print22, 'NUM = ',num
         gamma0=0.d0
         return
      end if
 22   format(10x,a6,i3)

c transfrom the energy and gamma values into quadruple precision
      print*, 'Input energy points'
      do i=1,num
         e_point(i)=e8(i)
         g_point(i)=g8(i)
         print*, e8(i),g8(i),i
      end do
      print*, ' '

c define the minimal energy and the maximal energies
      e_min=e_point(1)
      e_max=e_point(1)
      do i=2,num
         if( e_min.gt.e_point(i)) e_min=e_point(i)
         if( e_max.lt.e_point(i)) e_max=e_point(i)
      end do
c      print*, 'e_min =',e_min,'   e_max =',e_max

c define the energy of interest
      energy(1)=omega !epsiv+omega
c      if (energy(1)+epsiv.le.0.d0) then
c         print*, '***ERROR*** Stieltjes:'
c         print*, ' the photon energy is too small', energy(1),epsiv
c         return
c      end if

c check whether the defined energy belongs to the interval covered by E_POINT
      if (energy(1).lt.e_min.or.energy(1).gt.e_max) then
         gmean=0.d0
         do i=1,num
            gmean=gmean+g8(i)
         end do
         gmean=gmean/dfloat(num)
         print*, '***ERROR*** Stieltjes:'
         print*, ' the required energy is out of range'
         print33,energy(1),'is not in [',e_min,',',e_max,']'
         print*, ' the number of energy points is',num
         print*, ' the mean gamma in the covered range is',gmean
         print*, ' '
         gamma0=0.d0
         return
      end if
 33      format(10x,e16.10,2x,a11,e16.10,a1,e16.10,a1)

c shifting the energies to avoid the small denumerator problem
      do i=1,num
         e_point(i)=e_point(i)
c         print*, i,e_point(i)
      end do 

c initiate the recursive computation of the a,b coefficients and the orthogonal 
c polynomials according to (3.3.20-23) of Mueller-Plathe & Dierksen (1990)
       bcoef(0)=0.q0
       acoef(1)=0.q0
       do i=1,num
          bcoef(0)=bcoef(0)+g_point(i)
          acoef(1)=acoef(1)+g_point(i)/e_point(i)
       end do
       acoef(1)=acoef(1)/bcoef(0)

       do i=1,num
          qpol(0,i)=1.q0
          qpol(1,i)=1.q0/e_point(i)-acoef(1)
       end do

       bcoef(1)=0.q0
       acoef(2)=0.q0
       do i=1,num
          bcoef(1)=bcoef(1)+qpol(1,i)*g_point(i)/e_point(i)
          acoef(2)=acoef(2)+qpol(1,i)*g_point(i)/(e_point(i)**2)
       end do
       bcoef(1)= bcoef(1)/bcoef(0)
       acoef(2)=acoef(2)/(bcoef(0)*bcoef(1))-acoef(1)

c calculate the higher-order coefficients and polynomials recursively
c up to the (NUM-1)th order (total of NUM polynomials)
       asum=acoef(1)
       do i=3,num

          asum=asum+acoef(i-1)

          do j=1,num
             qpol(i-1,j)=(1.q0/e_point(j)-acoef(i-1))*qpol(i-2,j)
     &            -bcoef(i-2)*qpol(i-3,j)
          end do

          bprod=bcoef(0)
          do j=1,i-2
             bprod=bprod*bcoef(j)
          end do

          bcoef(i-1)=0.q0
          do j=1,num
             bcoef(i-1)=bcoef(i-1)+qpol(i-1,j)*g_point(j)
     &            /(e_point(j)**(i-1))
          end do
          bcoef(i-1)=bcoef(i-1)/bprod

          bprod=bprod*bcoef(i-1)

          acoef(i)=0.q0
          do j=1,num
             acoef(i)=acoef(i)+qpol(i-1,j)*g_point(j)/(e_point(j)**i)
          end do
          acoef(i)=acoef(i)/bprod-asum

c          print*, i,acoef(i),i-1,bcoef(i-1)

       end do

c calculate the NUM-th order polynomial just for the orthogonality check 
       do j=1,num
          qpol(num,j)=(1.q0/e_point(j)-acoef(num))*qpol(num-1,j)
     &         -bcoef(num-1)*qpol(num-2,j)
       end do

c check the orthogonality of the polynomials to define the maximal approximation order 
c if the orthogonality is preserved for all orders, MAXORD is set to the number of the 
c input points (NUM)
       maxord=num
       qnorm=bcoef(0)
       do i=1,num
          qnorm=0.q0
          qoverlap=0.q0
          do j=1,num 
             qnorm=qnorm+qpol(i,j)**2*g_point(j)
             qoverlap=qoverlap+qpol(i,j)*qpol(i-1,j)*g_point(j)
          end do
          if (qabs(qoverlap).lt.1.q-50) qoverlap=1.q-50
c          print*, i,qoverlap,qnorm
          if (qnorm/qabs(qoverlap).le.overmax) then
c MAXORD=I-1 is appropriate since the polynomial failing 
c the orthogonality check should not be used
             maxord=i-1
             go to 10
          end if
       end do

 10    continue

c look how many Stieltjes orders are available
       if (maxord.lt.5) then
          min=maxord
          max=maxord
          print*, '***WARNING*** Stieltjes:' 
          print*, ' only very low-order approximation is available'
          print*, ' MAXORD=',maxord
       else
          min=5
          max=maxord
       end if

c seek for convergence if enough Stieltjes orders are available
       converge=0
       if (maxord.ge.7) converge=1

c perform the gamma calculation using the successive approximations 
c N=5,...,MAXORD (if MAXORD > 5)
       do 20 imax=min,max

c fill the coefficients matrix
       do i=1,imax
          diag(i)=acoef(i)
       end do
       do i=2,imax
          offdiag(i)=-qsqrt(bcoef(i-1))
       end do

c diagonalize the coefficients matrix
c initialize the arrays
       do i=1,np
          do j=1,np
             abvec(i,j)=0.q0
          end do
          abvec(i,i)=1.q0
       end do
       call tql2(np,imax,diag,offdiag,abvec,ierr)
       if (ierr.ne.0) then
          print*, '***WARNING*** Stieltjes:'
          print*, ' the eigenvalue no. ',ierr,' failed to converge'
       end if
c fill the Stieltjes energy and gamma arrays
c note that the eigenvalues are inverse energies and are given in ascending order 
       do i=1,imax
          e_new(i)=1.q0/diag(imax+1-i)
          g_new(i)=bcoef(0)*abvec(1,imax+1-i)**2
c           print*, 'Die Punkte für die Diff'
c           print*, i,e_new(i),g_new(i)
       end do

c calculate the gamma's by simple numerical differentiation at the middle 
c point of each [E_NEW(I),E_NEW(I+1)] interval
       write(770+unit,*) imax
       do i=1,imax-1
          e0(i)=(e_new(i)+e_new(i+1))/2.d0
          gamma(i)=5.d-1*(g_new(i+1)+g_new(i))/(e_new(i+1)-e_new(i))
          write(770+unit,*) e0(i),gamma(i)
       end do
c       print*, ' ' 

c check whether the required energy is inside the interval covered by the E0 grid
      if (energy(1).lt.e0(1)) then 
         print*, '***WARNING*** Stieltjes:'
         print44, energy(1),' < ',e0(1),'(the first grid point)'
         print*,' large inaccuracy expected'
         gammas(1)=5.d-1*g_new(1)/e_new(1)
         print*, imax,gammas(1)
         gamma_h(imax)=gammas(1)
         go to 20
      end if
      if (energy(1).gt.e0(imax-1)) then 
         print*, '***WARNING*** Stieltjes:'
         print44, energy(1),' > ',e0(imax-1),'(the last grid point) '
         gammas(1)=0.d0
         print*, imax,gammas(1)
         gamma_h(imax)=gammas(1)
         go to 20
      end if
 44   format(10x,e16.10,a3,e16.10,x,a22)

c and interpolate the result using monotonicity-preserving piecewise 
c cubic Hermite interpolant 
       call e01bef (imax-1,e0,gamma,der,ifail)

c first for the photon energy of interest
       call e01bff (imax-1,e0,gamma,der,1,energy,gammas,ifail)
       gamma_h(imax)=gammas(1)
       print*, imax,gammas(1)
c print the gamma as a function of the Stieltjes order to a separate file
c       write(1000+77,*) imax,gammas(1)

c and then for the whole spectrum
       do i=1,nint-1
          energy_inter(i)=e0(1)+i*(e0(imax-1)-e0(1))/dfloat(nint)
       end do
       call e01bff (imax-1,e0,gamma,der,nint-1,energy_inter,
     &      gammas_inter,ifail)
       do i=1,nint-1
c          print*, energy_inter(i),gammas_inter(i),imax
c          write(2000+77,*) energy_inter(i),
c     &         gammas_inter(i),imax
       end do

c       write(2000+77,*) ' '

 20   continue

c      print*, ' '

c if not enough gamma values are available, the highest-order approximation is taken
      if (converge.eq.0) then
         gamma0=gamma_h(max)

c if enough gamma values are available to perform the convergence check, 
c the final gamma is taken to be an average of the three highest-order approximations
c unless convergence is detected; the converged result is preferred
      else if (converge.eq.1) then
         gamma0=(gamma_h(max-2)+gamma_h(max-1)+gamma_h(max))/3.d0
         conv_thresh=conv_thresh0 
         iconv=1

 40      continue

         do imax=max,min+2,-1

            dif12=dabs(gamma_h(imax)-gamma_h(imax-1))
            dif23=dabs(gamma_h(imax-1)-gamma_h(imax-2))
            dif13=dabs(gamma_h(imax)-gamma_h(imax-2))
            gmean=(gamma_h(imax-2)+gamma_h(imax-1)+gamma_h(imax))/3.d0
            difmean=(dif12+dif23+dif13)/3.d0

c prefer the convergence at higher orders over the one at lower orders
c to this end, the convergence criterion becomes more strict as we go to lower orders
            if (difmean/gmean.le.conv_thresh*conv_fac**(max-imax)) then
               print55, 'Convergence detected at',imax,'th order'
               print66, 'Convergence threshold =',conv_thresh*1.d2,'%'
               gamma0=gmean
               go to 30
            end if

         end do

c if we came here, the convergence was not found, so we increase the threshold
c maximal number of CONV_MAX convergence search loops is allowed
         iconv=iconv+1
         if (iconv.eq.conv_max) go to 30
         conv_thresh=conv_thresh*1.2d0
         go to 40 

      end if

 55   format(x,a23,x,i3,a8)
 66   format(x,a23,x,f4.1,a1)

 30   continue

      print*, 'The final result: Sigma =',gamma0
      print*, ' '
c      write(1000+77,*) ' '

      return
      end

