      subroutine bnd2td(nm,n,mb,a,d,e,e2,z)
      implicit none

      real*8 zero,half,one,two,dmin,dminrt,a,d,e,e2,z,g,b1,b2,s2,c2,
     .       f1,f2,u
      integer nm,n,mb,j,nm2,k,m1,n2,maxr,kr,r,r1,ugl,mr,j1,j2,l,i2,
     .        maxl,i1
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0,
     *           dmin=two**(-64),dminrt=two**(-32))
      dimension a(n,mb),d(n),e(n),e2(n),z(nm,n)

      do j=1,n
        d(j)=one
      enddo
      nm2=nm/2
      do j=1,n
        do k=1,nm
          z(k,j)=zero
        enddo
      enddo
      do j=1,nm2
        z(j,j)=one
        z(j+nm2,n-nm2+j)=one
      enddo

      m1=mb-1
      if (m1-1) 900,800,70

   70 n2=n-2

      do k=1,n2
        maxr=min(m1,n-k)
        do r1=2,maxr
          r=maxr+2-r1
          kr=k+r
          mr=mb-r
          g=a(kr,mr)
          a(kr-1,1)=a(kr-1,mr+1)
          ugl=k
          do j=kr,n,m1
            j1=j-1
            j2=j1-1
            if (g.eq.zero) exit
            b1=a(j1,1)/g
            b2=b1*d(j1)/d(j)
            s2=one/(one+b1*b2)
            if (s2.lt.half) then
              b1=g/a(j1,1)
              b2=b1*d(j)/d(j1)
              c2=one-s2
              d(j1)=c2*d(j1)
              d(j)=c2*d(j)
              f1=two*a(j,m1)
              f2=b1*a(j1,mb)
              a(j,m1)=-b2*(b1*a(j,m1)-a(j,mb))-f2+a(j,m1)
              a(j1,mb)=b2*(b2*a(j,mb)+f1)+a(j1,mb)
              a(j,mb)=b1*(f2-f1)+a(j,mb)
              do l=ugl,j2
                i2=mb-j+l
                u=a(j1,i2+1)+b2*a(j,i2)
                a(j,i2)=-b1*a(j1,i2+1)+a(j,i2)
                a(j1,i2+1)=u
              enddo
              ugl=j
              a(j1,1)=a(j1,1)+b2*g
              if (j.ne.n) then
                maxl=min(m1,n-j1)
                do l=2,maxl
                  i1=j1+l
                  i2=mb-l
                  u=a(i1,i2)+b2*a(i1,i2+1)
                  a(i1,i2+1)=-b1*a(i1,i2)+a(i1,i2+1)
                  a(i1,i2)=u
                enddo
                i1=j+m1
                if (i1.le.n) g=b2*a(i1,1)
              endif
              do l=1,nm
                u=z(l,j1)+b2*z(l,j)
                z(l,j)=-b1*z(l,j1)+z(l,j)
                z(l,j1)=u
              enddo
            else
              u=d(j1)
              d(j1)=s2*d(j)
              d(j)=s2*u
              f1=two*a(j,m1)
              f2=b1*a(j,mb)
              u=b1*(f2-f1)+a(j1,mb)
              a(j,m1)=b2*(b1*a(j,m1)-a(j1,mb))+f2-a(j,m1)
              a(j1,mb)=b2*(b2*a(j1,mb)+f1)+a(j,mb)
              a(j,mb)=u
              do l=ugl,j2
                i2=mb-j+l
                u=b2*a(j1,i2+1)+a(j,i2)
                a(j,i2)=-a(j1,i2+1)+b1*a(j,i2)
                a(j1,i2+1)=u
              enddo
              ugl=j
              a(j1,1)=b2*a(j1,1)+g
              if (j.ne.n) then
                maxl=min(m1,n-j1)
                do l=2,maxl
                  i1=j1+l
                  i2=mb-l
                  u=b2*a(i1,i2)+a(i1,i2+1)
                  a(i1,i2+1)=-a(i1,i2)+b1*a(i1,i2+1)
                  a(i1,i2)=u
                enddo
                i1=j+m1
                if (i1.le.n) then
                  g=a(i1,1)
                  a(i1,1)=b1*a(i1,1)
                endif
              endif
              do l=1,nm
                u=b2*z(l,j1)+z(l,j)
                z(l,j)=-z(l,j1)+b1*z(l,j)
                z(l,j1)=u
              enddo
            endif
          enddo
        enddo
        if (mod(k,64).eq.0) then
          do j=k,n
            if (d(j).ge.dmin) cycle
            maxl=max(1,mb+1-j)
            do l=maxl,m1
              a(j,l)=dminrt*a(j,l)
            enddo
            if (j.ne.n) then
              maxl=min(m1,n-j)
              do l=1,maxl
                i1=j+l
                i2=mb-l
                a(i1,i2)=dminrt*a(i1,i2)
              enddo
            endif
            do l=1,nm
              z(l,j)=dminrt*z(l,j)
            enddo
            a(j,mb)=dmin*a(j,mb)
            d(j)=d(j)/dmin
          enddo
        endif
      enddo

  800 do j=2,n
        e(j)=sqrt(d(j))
      enddo
      do j=2,n
        do k=1,nm
          z(k,j)=e(j)*z(k,j)
        enddo
      enddo
      u=one
      do j=2,n
        a(j,m1)=u*e(j)*a(j,m1)
        u=e(j)
        e2(j)=a(j,m1)**2
        a(j,mb)=d(j)*a(j,mb)
        d(j)=a(j,mb)
        e(j)=a(j,m1)
      enddo
      d(1)=a(1,mb)
      e(1)=zero
      e2(1)=zero
      return

  900 do j=1,n
        d(j)=a(j,mb)
        e(j)=zero
        e2(j)=zero
      enddo

      return
      end
