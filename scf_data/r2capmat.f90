
      subroutine CAPM(R_m,R_n,k_m,K_n,Integr)


      implicit none

      real *8 R_m(4),R_n(4),integr,data     
      integer *4 i,k_m(3),k_n(3),k
      data=0.0000001
     
      write(17,*)(R_m(i),i=1,3)
      write(17,*)(k_m(i),i=1,3)
      write(17,*)R_m(4)
      
      write(17,*)(R_n(i),i=1,3)
      write(17,*)(k_n(i),i=1,3)
      write(17,*)R_n(4)
      call mod_r2x(integr, R_m,R_n,k_m,k_n)
!      if(integr.gt.data)then
!      write(18,*)integr
!      endif
      
      return
      end
      

!C
      subroutine mod_r2x(result,x1,x2,k1,k2)
!C
!C-----------------------------------------------------------------------
!C The original version (r2x.f) of this program was written by Uwe Riss.
!C It calculated the matrix elements of r**2 (i.e. a parabolic CAP 
!C centered at the origin) with respect to GTOs. Such a CAP significantly
!C perturbs that region in space where the electrons still interact 
!C with each other. This modified version introduces a new parameter array
!C "core(3)" which refers to the distance from the symmetry center of the CAP
!C up to which the CAP is equal to 0 (along the x-, y- and z-axis, 
!C respectively). The location of the symmetry center of the CAP "cen(3)" is 
!C stated more explicitly.
!C                                  12.08.1998      Robin Santra
!C------------------------------------------------------------------------
      implicit real*8     (a-h,o-z)
!C
      real*8     x1(4), x2(4)
      real*8     g(0:8), h(0:8), f(0:20), q(0:16)
      real*8     ovl(3), dip(3)
      real*8     result
      real*8     core(3), cen(3)
      real*8     cl, cu
      integer    k1(3), k2(3)
      integer    i
       EXTERNAL dgamma
      common /boxsize/ core
      common /symcent/cen
!C
      alpha = x1(4)
      beta  = x2(4)

!C
!C-----------------------------------------------------------------------
!C     Calculation of the matrix elements
!C-----------------------------------------------------------------------
!C
!C     --------------------
!C     <I> - term (overlap)
!C     --------------------
!C
      do k=1,3
         m  = k1(k)
         n  = k2(k)
         nm = n + m
         ao = x1(k)
         bo = x2(k)
!C
      do i=0,m
         g(i) = 0.d0
      enddo
      g(m) = 1.d0
!C
      do i=0,n
         h(i) = 0.d0
      enddo
      h(n) = 1.d0
!C
      call gauss2( ao, alpha,  bo, beta, d0, delta, fakt )

!C d0 is the center and delta the exponent of the new gaussian 
!C formed by the two GTOs.
!C
      call poltrans( g(0), ao, d0, m )
      call poltrans( h(0), bo, d0, n )
      call polmulti( f(0), g(0), m, h(0), n )
!C
      temp = 0.d0
      do i=0,nm
         call intgauss( wert, delta, i )
         temp = temp + f(i) * wert
      enddo
      ovl(k) = fakt * temp
!C
      enddo
!C
!C     ------------
!C     <CAP> - term
!C     ------------
!C
      do k=1,3
         m  = k1(k)
         n  = k2(k)
         nm = n + m
         ao = x1(k)
         bo = x2(k)
!C
      do i=0,m
         g(i) = 0.d0
      enddo
      g(m) = 1.d0
!C
      do i=0,n
         h(i) = 0.d0
      enddo
      h(n) = 1.d0
!C
      call gauss2( ao, alpha,  bo, beta, d0, delta, fakt )
!C
      cl = cen(k) - d0 - core(k)
      cu = cen(k) - d0 + core(k)
!      write(*,*)cl,'cl'
!C The CAP is equal to 0 on the compact interval [cl,cu] (cl<=cu).
!C
      call poltrans( g(0), ao, d0, m )
      call poltrans( h(0), bo, d0, n )
      call polmulti( f(0), g(0), m, h(0), n )
!C
      temp = 0.d0
      do i=0,nm
         call intmod_r2( wert, delta, i, cl, cu )
         temp = temp + f(i) * wert
      enddo
      dip(k) = fakt * temp
!*     print*,'0>',m,n,ovl(k),dip(k)
!*     print*,'1>',f(0),f(1),f(2)
!C
      enddo
!C
!C-----------------------------------------------------------------------
!C
      result = dip(1) * ovl(2) * ovl(3)& 
            + dip(2) * ovl(3) * ovl(1) &
            + dip(3) * ovl(1) * ovl(2) 

!      result=ovl(1) * ovl(2)*OVl(3)
!C
!C-----------------------------------------------------------------------
!C
!C
!C     ----------
!C     <0> - Term
!C     ----------
!C
      alpha = 2.d0 * alpha
      beta  = 2.d0 * beta

!C
      do k=1,3
         m  = 2 * k1(k)
         n  = 2 * k2(k)

!write(111,*) alpha,beta
!write(*,*) alpha,beta
if(alpha.lt.0.00000000001)stop
if(beta.lt.0.00000000001)stop
         call intgauss( wert1, alpha, m )
         call intgauss( wert2, beta,  n )
!C
         ovl(k) = wert1 * wert2 
!C
      enddo
!C
!C-----------------------------------------------------------------------
!C
     
     
!      result= result/ sqrt( ovl(1) * ovl(2) * ovl(3) )

!C
!C-----------------------------------------------------------------------
!C
      return
!C
      end



!C
      subroutine intmod_r2(result,a,k,cl,cu)
!C
!C-----------------------------------------------------------
!C This subroutine calculates the integral of the function
!C x**k*f(x)*exp(-a*x**2) from -infinity to +infinity.
!C
!C        (x - cl)**2 if x < cl
!C f(x) =    0        if cl <= x <= cu 
!C        (x - cu)**2 if x > cu
!C
!C It exploites the fact that this integral can be expressed 
!C in terms of the incomplete GAMMA-function.    
!C                          12.08.1998 Robin Santra
!C-----------------------------------------------------------
      implicit none 
!C
      real*8   a, result
      real*8   cl, cu, dl, du, ql, qu
      real*8   par, sg, term
      integer  k, i, j, ex
      real*8 DGAMMA, inc_gamma, xi
       EXTERNAL dgamma     
      sg = (-1.d0)**k
      ql = a*cl*cl
      qu = a*cu*cu
!C
      result = 0.d0
      do i = 0, 2
         par = dble(k + i + 1)/2.d0
         ex = 2 - i
         dl = cl**ex
         du = (-cu)**ex
         j = k + i
         term =  (&
                DGAMMA(par)*(sg*dl + du)&
              - sg*dl*inc_gamma(par,ql)*xi(-cl,j)&
              -    du*inc_gamma(par,qu)*xi( cu,j)&
                ) / a**par
         if (i .eq. 1) then
            term = 2.d0*term
         end if
         result = result + term
      end do
!C      
      result = 0.5d0 * result
!C
      return 
!C
      end

      function inc_gamma(a,x)
!C-------------------------------------------------------
!C Returns the incomplete GAMMA function gamma(a,x)
!C which is defined by the integral of the function
!C exp(-t)*t**(a-1) over t from 0 to x.
!C This program is a modified version of a routine
!C provided in Press et al.: Numerical Recipes in FORTRAN,
!C !Cambridge University Press (1992), pages 211 - 213.
!C                                12.08.1998 Robin Santra
!C-------------------------------------------------------
      real*8 a, inc_gamma, x
      real*8 gammcf, gamser
      real*8 DGAMMA
       EXTERNAL dgamma
      if (x .lt. 0.d0 .or. a .le. 0.d0) write(*,*)'bad arguments in inc_gamma'
!C
      if (x .lt. a+1.d0) then
!C Use the series representation
         call gser(gamser,a,x)
         inc_gamma = gamser
      else
!C Use the continued fraction representation
         call gcf(gammcf,a,x)
!C and take its complement
         inc_gamma = DGAMMA(a) - gammcf
      end if
      return
      end
!C
!C
!C      
      subroutine gser(gamser,a,x)
!C      
      integer itmax
      real*8 a, gamser, x, eps
      parameter(itmax = 20000, eps = 1.d-20)
      integer n
      real*8 ap, del, sum
       EXTERNAL dgamma
      if (x .le. 0.d0) then
         if (x .lt. 0.d0) write(*,*)'x < 0 in gser'
         gamser = 0.d0
         return
      end if
!C
      ap = a
      sum = 1.d0/a
      del = sum
      do n = 1, itmax
         ap = ap + 1.d0
         del = del*x/ap
         sum = sum + del
         if (dabs(del) .lt. dabs(sum)*eps) goto 1
      end do
        write(*,*)'a too large, itmax too small in gser'
 1    gamser = sum*dexp(-x + a*dlog(x))
      return
      end
!C
!C
!C
      subroutine gcf(gammcf,a,x)
      integer itmax
      real*8 a, gammcf, x, eps, fpmin
      parameter(itmax = 20000, eps = 1.d-20,& 
               fpmin = 1.d-90)
!C itmax is the maximum allowed number of iterations;
!C eps is the relative accuracy; fpmin is a number near
!C the smallest representable floating point number.
      integer i
      real*8 an, b, c, d, del, h
       EXTERNAL dgamma
      b = x + 1.d0 -a
      c = 1.d0/fpmin
      d = 1.d0/b
      h = d
      do i = 1, itmax
         an = -i*(i - a)
         b = b + 2.d0
         d = an*d + b
         if (dabs(d) .lt. fpmin) d = fpmin
         c = b + an/c
         if (dabs(c) .lt. fpmin) c = fpmin
         d = 1.d0/d
         del = d*c
         h = h*del
         if (dabs(del - 1.d0) .lt. eps) goto 1
      end do
     write(*,*)  'a too large, itmax too small in gcf'
 1    gammcf = dexp(-x + a*dlog(x))*h
      return
      end

!c
!c This function is needed to determine a vital sign in CAP box 
!c routines "intmod_r2, intmod_r4, intbox".
!c
!c
      function xi(z,i)
      real*8 xi, z
      integer i
      if (z .ge. 0.d0) then
         xi = 1.d0
      else
         xi = -(-1.d0)**i
      end if
      return
      end







