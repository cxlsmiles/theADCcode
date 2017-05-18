      subroutine tddiag(nm,n,d,e,z,ierr)
      implicit real*8 (a-h,o-z)

      real*8 machep
      parameter (machep=2.22045d-16)
      parameter (zero=0.0d0,one=1.0d0,mxiter=30)
      dimension d(n),e(n),z(nm,n)

      ierr=0
      if (n.eq.1) return
      do i=2,n
        e(i-1)=e(i)
      enddo
      f=zero
      b=zero
      e(n)=zero
      do l=1,n
        j=0
        h=machep*(abs(d(l))+abs(e(l)))
        if (b.lt.h) b=h
        do m=l,n
          if (abs(e(m)).le.b) exit
        enddo
        do while (abs(e(l)).gt.b.and.j.lt.mxiter)
          j=j+1
          l1=l+1
          g=d(l)
          p=(d(l1)-g)/(e(l)*2)
          r=sqrt(p**2+one)
          d(l)=e(l)/(p+sign(r,p))
          h=g-d(l)
          do i=l1,n
            d(i)=d(i)-h
          enddo
          f=f+h
          p=d(m)
          c=one
          s=zero
          mml=m-l
          do ii=1,mml
            i=m-ii
            g=c*e(i)
            h=c*p
            if (abs(p).ge.abs(e(i))) then
              c=e(i)/p
              r=sqrt(c**2+one)
              e(i+1)=s*p*r
              s=c/r
              c=one/r
            else
              c=p/e(i)
              r=sqrt(c**2+one)
              e(i+1)=s*e(i)*r
              s=one/r
              c=c*s
            endif
            p=c*d(i)-s*g
            d(i+1)=h+s*(c*g+s*d(i))
            do k=1,nm
              h=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*h
              z(k,i)=c*z(k,i)-s*h
            enddo
          enddo
          e(l)=s*p
          d(l)=c*p
        enddo
        if (j.eq.mxiter) goto 1000
        d(l)=d(l)+f
      enddo
      do ii=2,n
        i=ii-1
        k=i
        p=d(i)
        do j=ii,n
          if (d(j).ge.p) cycle
          k=j
          p=d(j)
        enddo
        if (k.eq.i) cycle
        d(k)=d(i)
        d(i)=p
        do j=1,nm
          p=z(j,i)
          z(j,i)=z(j,k)
          z(j,k)=p
        enddo
      enddo
      return
 1000 ierr=l
      return
      end
