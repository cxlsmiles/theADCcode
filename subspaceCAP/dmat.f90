module dmat
use CAP_MO
implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function vabcd(a,b,c,d)
   implicit none

    double precision :: vabcd

    integer, intent(in) :: a, b, c, d

    double precision :: vpqrs
    external vpqrs

! nico 16.06.2011
! minus comes from different definitions of higher order pertub. corrections of the wf. see Yasen.
      vabcd = -vpqrs(a,c,b,d)

   end function vabcd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function tmat(it,i,j,ip,jp)
   implicit none

    double precision :: tmat
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a

      tmat = 0d0

    if(it==1) then

      if(i==ip) then
       do a = nocc+1, nMO
        tmat = tmat + DMO(a,i)*mat_CAP_MO(ip,a)
       enddo 
        tmat = -4d0*tmat
      endif

    elseif(it==2) then

      if(i==ip) then
       do a = nocc+1, nMO
          tmat = tmat + DMO(a,j)*mat_CAP_MO(ip,a) + DMO(a,ip)*mat_CAP_MO(j,a)
       enddo 
      endif

      if(j==ip) then
       do a = nocc+1, nMO
          tmat = tmat + DMO(a,i)*mat_CAP_MO(ip,a)  + DMO(a,ip)*mat_CAP_MO(i,a)
       enddo 
      endif

        tmat = -2d0*tmat/dsqrt(2d0)

    elseif(it==3) then

      if(i==ip) then
       do a = nocc+1, nMO
          tmat = tmat + DMO(a,j)*mat_CAP_MO(jp,a)  &
        + DMO(a,jp)*mat_CAP_MO(j,a)
       enddo 
      endif

      if(j==jp) then
       do a = nocc+1, nMO
          tmat = tmat + DMO(a,i)*mat_CAP_MO(ip,a) &
        + DMO(a,ip)*mat_CAP_MO(i,a)
       enddo 
      endif

      if(i==jp) then
       do a = nocc+1, nMO
          tmat = tmat + DMO(a,j)*mat_CAP_MO(ip,a) &
        + DMO(a,ip)*mat_CAP_MO(j,a)
       enddo 
      endif

      if(j==ip) then
       do a = nocc+1, nMO
          tmat = tmat + DMO(a,i)*mat_CAP_MO(jp,a) &
        + DMO(a,jp)*mat_CAP_MO(i,a)
       enddo 
      endif

          tmat = -1d0*tmat

    endif

   end function  tmat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d21(it,i,j,ip,jp)
   implicit none

    double precision :: d21
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, ap, b, bp, k
    double precision :: t

!    !double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d21 = 0d0
     if(i/=ip) return

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do ap = nocc + 1, nMO
           do k = 1, nocc
                d21 = d21 + (2d0*vabcd(ap,b,k,jp)*vabcd(a,b,k,j) &
                            +2d0*vabcd(ap,b,jp,k)*vabcd(a,b,j,k) &
                            -1d0*vabcd(ap,b,jp,k)*vabcd(a,b,k,j) &
                            -1d0*vabcd(ap,b,k,jp)*vabcd(a,b,j,k)) &
                            *mat_CAP_MO(ap,a)/                        &
                            ((EMO(ap)+EMO(b)-EMO(k)-EMO(jp))*(EMO(a)+EMO(b)-EMO(k)-EMO(j)))
           enddo
         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
        d21 = - 0.5 * d21            
     elseif(i==j .or. ip==jp) then
        d21 = - 0.5 * d21/dsqrt(2d0)            
     else
        d21 = - 0.5 * d21
     endif

   end function d21


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d22(it,i,j,ip,jp)
   implicit none

    double precision :: d22
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k, kp
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d22 = 0d0
    if(i/=ip) return


     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc
           do kp = 1, nocc

                d22 = d22 + (2d0*vabcd(a,b,kp,k)*vabcd(a,b,j,k) &
                            +2d0*vabcd(a,b,k,kp)*vabcd(a,b,k,j) &
                            -1d0*vabcd(a,b,kp,k)*vabcd(a,b,k,j) &
                            -1d0*vabcd(a,b,k,kp)*vabcd(a,b,j,k)) &
                            *mat_CAP_MO(jp,kp)/                        &
                            ((EMO(a)+EMO(b)-EMO(k)-EMO(kp))*(EMO(a)+EMO(b)-EMO(k)-EMO(j))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
        d22 = 0.25d0*d22            
     elseif(i==j .or. ip==jp) then
        d22 = 0.25d0*d22/dsqrt(2d0)            
     else
        d22 = 0.25d0*d22
     endif

   end function d22

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d23(it,i,j,ip,jp)
   implicit none

    double precision :: d23
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k, kp
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d23 = 0d0
    if(i/=ip) return

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc
           do kp = 1, nocc

                d23 = d23 + (2d0*vabcd(a,b,jp,kp)*vabcd(a,b,j,k) &
                            +2d0*vabcd(a,b,kp,jp)*vabcd(a,b,k,j) &
                            -1d0*vabcd(a,b,jp,kp)*vabcd(a,b,k,j) &
                            -1d0*vabcd(a,b,kp,jp)*vabcd(a,b,j,k)) &
                            *mat_CAP_MO(k,kp)/                        &
                            ((EMO(a)+EMO(b)-EMO(kp)-EMO(jp))*(EMO(a)+EMO(b)-EMO(k)-EMO(j))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
     d23 = 0.25d0*d23            
     elseif(i==j .or. ip==jp) then
     d23 = 0.25d0*d23/dsqrt(2d0)            
     else
     d23 = 0.25d0*d23
     endif


   end function d23

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d24(it,i,j,ip,jp)
   implicit none

    double precision :: d24
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, ap, k
    double precision :: t


!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d24 = 0d0
     if(j/=jp) return

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do ap = nocc + 1, nMO
           do k = 1, nocc

                d24 = d24 + (2d0*vabcd(ap,b,k,ip)*vabcd(a,b,k,i) &
                            +2d0*vabcd(ap,b,ip,k)*vabcd(a,b,i,k) &
                            -1d0*vabcd(ap,b,ip,k)*vabcd(a,b,k,i) &
                            -1d0*vabcd(ap,b,k,ip)*vabcd(a,b,i,k)) &
                            *mat_CAP_MO(ap,a)/                        &
                            ((EMO(ap)+EMO(b)-EMO(k)-EMO(ip))*(EMO(a)+EMO(b)-EMO(k)-EMO(i))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
     d24 = - 0.5 * d24            
     elseif(i==j .or. ip==jp) then
     d24 = - 0.5 * d24/dsqrt(2d0)            
     else
     d24 = - 0.5 * d24
     endif

   end function d24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d25(it,i,j,ip,jp)
   implicit none

    double precision :: d25
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k, kp
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d25 = 0d0
     if(j/=jp) return

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc
           do kp = 1, nocc

                d25 = d25 + (2d0*vabcd(a,b,k,kp)*vabcd(a,b,k,i) &
                            +2d0*vabcd(a,b,kp,k)*vabcd(a,b,i,k) &
                            -1d0*vabcd(a,b,k,kp)*vabcd(a,b,i,k) &
                            -1d0*vabcd(a,b,kp,k)*vabcd(a,b,k,i)) &
                            *mat_CAP_MO(ip,kp)/                        &
                            ((EMO(a)+EMO(b)-EMO(k)-EMO(kp))*(EMO(a)+EMO(b)-EMO(k)-EMO(i))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
     d25 = 0.25d0*d25            
     elseif(i==j .or. ip==jp) then
     d25 = 0.25d0*d25/dsqrt(2d0)            
     else
     d25 = 0.25d0*d25
     endif

   end function d25

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d26(it,i,j,ip,jp)
   implicit none

    double precision :: d26
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k, kp
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d26 = 0d0
     if(j/=jp) return

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc
           do kp = 1, nocc

                d26 = d26 + (2d0*vabcd(a,b,kp,ip)*vabcd(a,b,k,i) &
                            +2d0*vabcd(a,b,ip,kp)*vabcd(a,b,i,k) &
                            -1d0*vabcd(a,b,kp,ip)*vabcd(a,b,i,k) &
                            -1d0*vabcd(a,b,ip,kp)*vabcd(a,b,k,i)) &
                            *mat_CAP_MO(k,kp)/                        &
                            ((EMO(a)+EMO(b)-EMO(kp)-EMO(ip))*(EMO(a)+EMO(b)-EMO(k)-EMO(i))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
       d26 = 0.25d0*d26            
     elseif(i==j .or. ip==jp) then
       d26 = 0.25d0*d26/dsqrt(2d0)            
     else
       d26 = 0.25d0*d26
     endif

   end function d26

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d27(it,i,j,ip,jp)
   implicit none

    double precision :: d27
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, ap, k
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d27 = 0d0
     if(i/=jp) return
     if(i==j .and. ip==jp) then
       d27 = 0d0            
       return
     endif

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do ap = nocc + 1, nMO
           do k = 1, nocc

                d27 = d27 + (2d0*vabcd(ap,b,k,ip)*vabcd(a,b,k,j) &
                            +2d0*vabcd(ap,b,ip,k)*vabcd(a,b,j,k) &
                            -1d0*vabcd(ap,b,k,ip)*vabcd(a,b,j,k) &
                            -1d0*vabcd(ap,b,ip,k)*vabcd(a,b,k,j)) &
                            *mat_CAP_MO(ap,a)/                        &
                            ((EMO(ap)+EMO(b)-EMO(k)-EMO(ip))*(EMO(a)+EMO(b)-EMO(k)-EMO(j))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .or. ip==jp) then
       d27 = - 0.5 * d27/dsqrt(2d0)            
     else
       d27 = - 0.5 * d27
     endif

   end function d27

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d28(it,i,j,ip,jp)
   implicit none

    double precision :: d28
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k, kp
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d28 = 0d0
     if(i/=jp) return
     if(i==j .and. ip==jp) then
       d28 = 0d0
       return
     endif

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc
           do kp = 1, nocc

                d28 = d28 + (2d0*vabcd(a,b,k,kp)*vabcd(a,b,k,j) &
                            +2d0*vabcd(a,b,kp,k)*vabcd(a,b,j,k) &
                            -1d0*vabcd(a,b,k,kp)*vabcd(a,b,j,k) &
                            -1d0*vabcd(a,b,kp,k)*vabcd(a,b,k,j)) &
                            *mat_CAP_MO(ip,kp)/                        &
                            ((EMO(a)+EMO(b)-EMO(k)-EMO(kp))*(EMO(a)+EMO(b)-EMO(k)-EMO(j))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .or. ip==jp) then
     d28 = 0.25d0*d28/dsqrt(2d0)            
     else
     d28 = 0.25d0*d28
     endif

   end function d28

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d29(it,i,j,ip,jp)
   implicit none

    double precision :: d29
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k, kp
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d29 = 0d0
     if(i/=jp) return

     if(i==j .and. ip==jp) then
       d29 = 0d0
       return
     endif

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc
           do kp = 1, nocc

                d29 = d29 + (2d0*vabcd(a,b,kp,ip)*vabcd(a,b,k,j) &
                            +2d0*vabcd(a,b,ip,kp)*vabcd(a,b,j,k) &
                            -1d0*vabcd(a,b,kp,ip)*vabcd(a,b,j,k) &
                            -1d0*vabcd(a,b,ip,kp)*vabcd(a,b,k,j)) &
                            *mat_CAP_MO(k,kp)/                        &
                            ((EMO(a)+EMO(b)-EMO(kp)-EMO(ip))*(EMO(a)+EMO(b)-EMO(k)-EMO(j))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .or. ip==jp) then
     d29 = 0.25d0*d29/dsqrt(2d0)            
     else
     d29 = 0.25d0*d29
     endif

   end function d29

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d210(it,i,j,ip,jp)
   implicit none

    double precision :: d210
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, ap, k
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d210 = 0d0
     if(j/=ip) return
     if(i==j .and. ip==jp) then
       d210 = 0d0
       return
     endif

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do ap = nocc + 1, nMO
           do k = 1, nocc

                d210 = d210 + (2d0*vabcd(ap,b,k,jp)*vabcd(a,b,k,i) &
                            +2d0*vabcd(ap,b,jp,k)*vabcd(a,b,i,k) &
                            -1d0*vabcd(ap,b,k,jp)*vabcd(a,b,i,k) &
                            -1d0*vabcd(ap,b,jp,k)*vabcd(a,b,k,i)) &
                            *mat_CAP_MO(ap,a)/                        &
                            ((EMO(ap)+EMO(b)-EMO(k)-EMO(jp))*(EMO(a)+EMO(b)-EMO(k)-EMO(i))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .or. ip==jp) then
       d210 = - 0.5 * d210/dsqrt(2d0)            
     else
       d210 = - 0.5 * d210
     endif

   end function d210

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d211(it,i,j,ip,jp)
   implicit none

    double precision :: d211
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k, kp
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d211 = 0d0
     if(j/=ip) return
     if(i==j .and. ip==jp) then
       d211 = 0d0
       return
     endif

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc
           do kp = 1, nocc

                d211 = d211 + (2d0*vabcd(a,b,kp,k)*vabcd(a,b,i,k) &
                            +2d0*vabcd(a,b,k,kp)*vabcd(a,b,k,i) &
                            -1d0*vabcd(a,b,kp,k)*vabcd(a,b,k,i) &
                            -1d0*vabcd(a,b,k,kp)*vabcd(a,b,i,k)) &
                            *mat_CAP_MO(jp,kp)/                        &
                            ((EMO(a)+EMO(b)-EMO(k)-EMO(kp))*(EMO(a)+EMO(b)-EMO(k)-EMO(i))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .or. ip==jp) then
     d211 = 0.25d0*d211/dsqrt(2d0)            
     else
     d211 = 0.25d0*d211
     endif

   end function d211

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d212(it,i,j,ip,jp)
   implicit none

    double precision :: d212
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k, kp
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d212 = 0d0
     if(j/=ip) return
     if(i==j .and. ip==jp) then
       d212 = 0d0
       return
     endif

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc
           do kp = 1, nocc

                d212 = d212 + (2d0*vabcd(a,b,kp,jp)*vabcd(a,b,k,i) &
                            +2d0*vabcd(a,b,jp,kp)*vabcd(a,b,i,k) &
                            -1d0*vabcd(a,b,kp,jp)*vabcd(a,b,i,k) &
                            -1d0*vabcd(a,b,jp,kp)*vabcd(a,b,k,i)) &
                            *mat_CAP_MO(k,kp)/                        &
                            ((EMO(a)+EMO(b)-EMO(kp)-EMO(jp))*(EMO(a)+EMO(b)-EMO(k)-EMO(i))) 
           enddo
         enddo
       enddo
     enddo

     if(i==j .or. ip==jp) then
     d212 = 0.25d0*d212/dsqrt(2d0)            
     else
     d212 = 0.25d0*d212
     endif

   end function d212

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d213(it,i,j,ip,jp)
   implicit none

    double precision :: d213
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, ap
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d213 = 0d0

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do ap = nocc + 1, nMO

                d213 = d213 + (vabcd(ap,b,ip,jp)*vabcd(a,b,i,j)+vabcd(ap,b,jp,ip)*vabcd(a,b,j,i) &
                              +vabcd(ap,b,ip,jp)*vabcd(a,b,j,i)+vabcd(ap,b,jp,ip)*vabcd(a,b,i,j)) &
                               *mat_CAP_MO(ap,a)/((EMO(ap)+EMO(b)-EMO(ip)-EMO(jp))*(EMO(a)+EMO(b)-EMO(i)-EMO(j)))

         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
        d213 = 0.25 * d213
     elseif(i==j .or. ip==jp) then
        d213 = 0.5 * d213/dsqrt(2d0)
     else         
        d213 = 0.5 * d213
     endif

   end function d213

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d214(it,i,j,ip,jp)
   implicit none

    double precision :: d214
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d214 = 0d0

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc

                d214 = d214 + (vabcd(a,b,k,ip)*vabcd(a,b,j,i)+vabcd(a,b,ip,k)*vabcd(a,b,i,j) &
                              +vabcd(a,b,k,ip)*vabcd(a,b,i,j)+vabcd(a,b,ip,k)*vabcd(a,b,j,i)) &
                               *mat_CAP_MO(jp,k)/((EMO(a)+EMO(b)-EMO(ip)-EMO(k))*(EMO(a)+EMO(b)-EMO(i)-EMO(j)))

         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
        d214 = -0.125d0*d214
     elseif(i==j .or. ip==jp) then
        d214 = -0.25d0*d214/dsqrt(2d0)            
     else         
        d214 = -0.25d0*d214
     endif

   end function d214

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function d215(it,i,j,ip,jp)
   implicit none

    double precision :: d215
    integer, intent(in) :: it, i, j, ip, jp
    integer :: a, b, k
    double precision :: t

!    double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     d215 = 0d0

     do a = nocc + 1, nMO
       do b = nocc + 1, nMO
         do k = 1, nocc

                d215 = d215 + (vabcd(a,b,jp,k)*vabcd(a,b,j,i)+vabcd(a,b,k,jp)*vabcd(a,b,i,j)  &
                              +vabcd(a,b,jp,k)*vabcd(a,b,i,j)+vabcd(a,b,k,jp)*vabcd(a,b,j,i)) &
                               *mat_CAP_MO(ip,k)/((EMO(a)+EMO(b)-EMO(k)-EMO(jp))*(EMO(a)+EMO(b)-EMO(i)-EMO(j)))

         enddo
       enddo
     enddo

     if(i==j .and. ip==jp) then
        d215 = -0.125d0*d215
     elseif(i==j .or. ip==jp) then
        d215 = -0.25d0*d215/dsqrt(2d0)            
     else         
        d215 = -0.25d0*d215
     endif


   end function d215

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function t21(i,j,ap,ip,jp,kp)
   implicit none

    double precision :: t21
    integer, intent(in) :: i, j, ap, ip, jp, kp
    integer :: b, k

    !double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     t21 = 0d0

    if(i==jp) then

     do b = nocc+1, nMO
       do k = 1, nocc
         t21 = t21 + (4d0*vabcd(ap,b,ip,k)-2d0*vabcd(ap,b,k,ip))*mat_CAP_MO(b,k)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(k))                              
       enddo
     enddo

     do b = nocc+1, nMO
         t21 = t21 - (4d0*vabcd(ap,b,ip,jp)-2d0*vabcd(ap,b,jp,ip))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(jp))
     enddo
   endif

    if(i==ip) then
     do b = nocc+1, nMO
         t21 = t21 + 2d0*vabcd(ap,b,jp,jp)*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(jp))
     enddo
   endif

         t21 = -t21/dsqrt(2d0)
!         write(*,*)'NICO',ap,i,jp,t21

   end function t21

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function t22(i,j,ap,ip,jp,kp)
   implicit none

    double precision :: t22
    integer, intent(in) ::  i, j, ap, ip, jp, kp
    integer :: b, k

    !double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     t22 = 0d0

    if(i==ip .and. j==jp) then
     do b = nocc+1, nMO
       do k = 1, nocc
         t22 = t22 + (2d0*vabcd(ap,b,jp,k)-vabcd(ap,b,k,jp))*mat_CAP_MO(b,k)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(k))
       enddo
     enddo
    endif

    if(i==jp .and. j==ip) then
     do b = nocc+1, nMO
       do k = 1, nocc
         t22 = t22 + (2d0*vabcd(ap,b,jp,k)-vabcd(ap,b,k,jp))*mat_CAP_MO(b,k)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(k))
       enddo
     enddo
    endif

    if(i==ip) then
     do b = nocc+1, nMO
         t22 = t22 - vabcd(ap,b,jp,jp)*mat_CAP_MO(b,j)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(jp))
     enddo
    endif

    if(j==ip) then
     do b = nocc+1, nMO
         t22 = t22 - vabcd(ap,b,jp,jp)*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(jp))
     enddo
    endif

    if(i==jp) then
     do b = nocc+1, nMO
         t22 = t22 + (2d0*vabcd(ap,b,ip,jp)-vabcd(ap,b,jp,ip))*mat_CAP_MO(b,j)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(ip))
     enddo
    endif

    if(j==jp) then
     do b = nocc+1, nMO
         t22 = t22 + (2d0*vabcd(ap,b,ip,jp)-vabcd(ap,b,jp,ip))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(ip))
     enddo
    endif

   end function t22

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function t23(i,j,ap,ip,jp,kp)
   implicit none

    double precision :: t23
    integer, intent(in) ::  i, j, ap, ip, jp, kp
    integer :: b, k

    !double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     t23 = 0d0

    if(i==ip) then
     do b = nocc+1, nMO
         t23 = t23 + (2d0*vabcd(ap,b,jp,kp)-vabcd(ap,b,kp,jp))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(jp))
     enddo
    endif

    if(i==jp) then
     do b = nocc+1, nMO
         t23 = t23 - (vabcd(ap,b,kp,ip)+vabcd(ap,b,ip,kp))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(ip))
     enddo
    endif

    if(i==kp) then
     do b = nocc+1, nMO
         t23 = t23 + (2d0*vabcd(ap,b,jp,ip)-vabcd(ap,b,ip,jp))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(jp))
     enddo
    endif

   end function t23

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function t24(i,j,ap,ip,jp,kp)
   implicit none

    double precision :: t24
    integer, intent(in) :: i, j, ap, ip, jp, kp
    integer :: b, k

    double precision :: temp

    !double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     t24 = 0d0
     temp = 0d0

    if(i==ip .and. j==jp) then
     do b = nocc+1, nMO
       do k = 1, nocc
         t24 = t24 + (2d0*vabcd(ap,b,kp,k)-vabcd(ap,b,k,kp))*mat_CAP_MO(b,k)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(k))
       enddo
     enddo
    endif

        temp=t24

    if(i==jp .and. j==kp) then
     do b = nocc+1, nMO
       do k = 1, nocc
         t24 = t24 + (2d0*vabcd(ap,b,ip,k)-vabcd(ap,b,k,ip))*mat_CAP_MO(b,k)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(k))
       enddo
     enddo
    endif
        temp=t24

    if(i==ip .and. j==kp) then
     do b = nocc+1, nMO
       do k = 1, nocc
         t24 = t24 - 2d0*(2d0*vabcd(ap,b,jp,k)-vabcd(ap,b,k,jp))*mat_CAP_MO(b,k)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(k))
       enddo
     enddo
    endif

        temp=t24

    if(i==ip) then
     do b = nocc+1, nMO
         t24 = t24 + (2d0*vabcd(ap,b,jp,kp)-vabcd(ap,b,kp,jp))*mat_CAP_MO(b,j)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(jp))
     enddo
    endif

        temp=t24

    if(j==ip) then
     do b = nocc+1, nMO
         t24 = t24 + (2d0*vabcd(ap,b,jp,kp)-vabcd(ap,b,kp,jp))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(jp))
     enddo
    endif

        temp=t24

    if(i==jp) then
     do b = nocc+1, nMO
         t24 = t24 - (vabcd(ap,b,kp,ip)+vabcd(ap,b,ip,kp))*mat_CAP_MO(b,j)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(ip))
     enddo
    endif

        temp=t24

    if(j==jp) then
     do b = nocc+1, nMO
         t24 = t24 - (vabcd(ap,b,kp,ip)+vabcd(ap,b,ip,kp))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(ip))
     enddo
    endif

        temp=t24

    if(i==kp) then
     do b = nocc+1, nMO
         t24 = t24 + (2d0*vabcd(ap,b,jp,ip)-vabcd(ap,b,ip,jp))*mat_CAP_MO(b,j)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(jp))
     enddo
    endif

        temp=t24

    if(j==kp) then
     do b = nocc+1, nMO
         t24 = t24 + (2d0*vabcd(ap,b,jp,ip)-vabcd(ap,b,ip,jp))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(jp))
     enddo
    endif

        temp=t24

         t24 = t24/dsqrt(2d0)

   end function t24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function t25(i,j,ap,ip,jp,kp)
   implicit none

    double precision :: t25
    integer, intent(in) :: i, j, ap, ip, jp, kp
    integer :: b, k

    !double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     t25 = 0d0

    if(i==ip) then
     do b = nocc+1, nMO
         t25 = t25 + vabcd(ap,b,kp,jp)*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(jp))
     enddo
    endif

    if(i==jp) then
     do b = nocc+1, nMO
         t25 = t25 + (vabcd(ap,b,kp,ip)-vabcd(ap,b,ip,kp))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(ip))
     enddo
    endif

    if(i==kp) then
     do b = nocc+1, nMO
         t25 = t25 - vabcd(ap,b,ip,jp)*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(jp)-EMO(ip))
     enddo
    endif

         t25 = dsqrt(3d0) * t25

   end function t25

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function t26(i,j,ap,ip,jp,kp)
   implicit none

    double precision :: t26
    integer, intent(in) :: i, j, ap, ip, jp, kp
    integer :: b, k

    !double precision :: vabcd
    double precision :: vpqrs
    external vpqrs

     t26 = 0d0

    if(i==ip .and. j==jp) then
     do b = nocc+1, nMO
       do k = 1, nocc
         t26 = t26 - (2d0*vabcd(ap,b,kp,k)-vabcd(ap,b,k,kp))*mat_CAP_MO(b,k)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(k))
       enddo
     enddo
    endif

    if(i==jp .and. j==kp) then
     do b = nocc+1, nMO
       do k = 1, nocc
         t26 = t26 + (2d0*vabcd(ap,b,ip,k)-vabcd(ap,b,k,ip))*mat_CAP_MO(b,k)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(k))
       enddo
     enddo
    endif

    if(i==ip) then
     do b = nocc+1, nMO
         t26 = t26 + vabcd(ap,b,kp,jp)*mat_CAP_MO(b,j)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(jp))
     enddo
    endif

    if(j==ip) then
     do b = nocc+1, nMO
         t26 = t26 + vabcd(ap,b,kp,jp)*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(jp))
     enddo
    endif


    if(i==jp) then
     do b = nocc+1, nMO
         t26 = t26 + (vabcd(ap,b,kp,ip)-vabcd(ap,b,ip,kp))*mat_CAP_MO(b,j)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(ip))
     enddo
    endif

    if(j==jp) then
     do b = nocc+1, nMO
         t26 = t26 + (vabcd(ap,b,kp,ip)-vabcd(ap,b,ip,kp))*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(kp)-EMO(ip))
     enddo
    endif

    if(i==kp) then
     do b = nocc+1, nMO
         t26 = t26 - vabcd(ap,b,ip,jp)*mat_CAP_MO(b,j)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(jp))
     enddo
    endif

    if(j==kp) then
     do b = nocc+1, nMO
         t26 = t26 - vabcd(ap,b,ip,jp)*mat_CAP_MO(b,i)/ &
                     (EMO(ap)+EMO(b)-EMO(ip)-EMO(jp))
     enddo
    endif

         t26 = dsqrt(1.5d0) * t26

   end function t26

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine test_me(nmo,dens)
    integer, intent(in) :: nmo
    double precision, intent(in):: dens(nmo,nmo)
    
    integer :: i,j

    double precision :: vpqrs
    external vpqrs

    write(*,*) 'this is a test'
    write(*,*) 'print one integral:', vpqrs(1,1,1,1)
    write(*,*) 'print the density matrix'

    do i=1,nmo
       do j=1,nmo
          write(*,*) i, j, dens(i,j)
       enddo
    enddo
  end subroutine test_me
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module dmat
