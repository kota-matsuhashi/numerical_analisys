C----------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(if=21,jf=21,kf=21)
      common T(if,jf,kf)
      common dx,dy,dz
      dimension al(if,jf,kf)
C
      data omega/1.5/
      data ia,ib/8,14/
      data ja,jb/8,14/
      data ka,kb/8,14/
      data ck1,ck2/10.0,1.0/
C
      dt=0.001
      dx=0.1
      dy=0.1
      dz=0.1
      dx2=dx*dx
      dy2=dy*dy
      dz2=dz*dz
      dlx=dt/dx2
      dly=dt/dy2
      dlz=dt/dz2
C
      do k=1,kf
      do j=1,jf
      do i=1,if
        T(i,j,k)=50.0
        al(i,j,k)=50.0
      enddo
      enddo
      enddo
C
      do k=1,kf
      do j=1,jf
        T(1,j,k)=100.0
        T(if,j,k)=50.0
      enddo
      enddo
C
      iter=0
      tsum=0.0
   99 iter=iter+1
      tsum=tsum+dt
      do i=2,if-1
      do j=1,jf
      do k=1,kf
c
      im1=i-1
      if (i.eq.1) im1=2
      ip1=i+1
      if (i.eq.if) ip1=if-1
      ckxl=ck1
      ckxr=ck1
      if ((j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (i.gt.ia).and.(i.le.ib)) ckxl=ck2
      if ((j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (i.ge.ia).and.(i.lt.ib)) ckxr=ck2
c
      jm1=j-1
      if (j.eq.1) jm1=2
      jp1=j+1
      if (j.eq.jf) jp1=jf-1
      ckyl=ck1
      ckyr=ck1
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (j.gt.ja).and.(j.le.jb)) ckyl=ck2
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (j.ge.ja).and.(j.lt.jb)) ckyr=ck2
c
      km1=k-1
      if (k.eq.1) km1=2
      kp1=k+1
      if (k.eq.kf) kp1=kf-1
      ckzl=ck1
      ckzr=ck1
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (j.ge.ja).and.(j.le.jb).and.
     &    (k.gt.ka).and.(k.le.kb)) ckzl=ck2
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.lt.kb)) ckzr=ck2
C
      r=(ckxl*(T(im1,j,k)-T(i,j,k))
     &  +ckxr*(T(ip1,j,k)-T(i,j,k)))*dlx*0.5
     & +(ckyl*(T(i,jm1,k)-T(i,j,k))
     &  +ckyr*(T(i,jp1,k)-T(i,j,k)))*dly*0.5
     & +(ckzl*(T(i,j,km1)-T(i,j,k))
     &  +ckzr*(T(i,j,kp1)-T(i,j,k)))*dlz*0.5
     &  +T(i,j,k)
      al(i,j,k)=r
      enddo
      enddo
      enddo
c
      na=0
  100 na=na+1
      rmax=0.0
C
      do ii=2,if-1
      i=ii
      if (mod(na,2).eq.0) i=if+1-i
      do jj=1,jf
      j=jj
      if (mod(na,2).eq.0) j=jf+1-j
      do kk=1,kf
      k=kk
      if (mod(na,2).eq.0) k=kf+1-k
C
      dl=1.0+(dlx*(ckxr+ckxl)
     &       +dly*(ckyr+ckyl)
     &       +dlz*(ckzr+ckzl))*0.5
      r=(ckxl*T(im1,j,k)+ckxr*T(ip1,j,k))*dlx/(dl*2.0)
     & +(ckyl*T(i,jm1,k)+ckyr*T(i,jp1,k))*dly/(dl*2.0)
     & +(ckzl*T(i,j,km1)+ckzr*T(i,j,kp1))*dlz/(dl*2.0)
     &  +al(i,j,k)/dl - T(i,j,k)
      T(i,j,k)=T(i,j,k) + omega*r
      rmax=amax1(rmax,abs(r))
      enddo
      enddo
      enddo
c
      if (na.eq.1) rmax1=rmax
      rmax=rmax/rmax1
      write(6,20) na, rmax
      write(11,20) na,rmax
   20 format(5x,'na=',i3,5x,'rmax=',f10.7)
c
      if (na.lt.10.and.rmax.gt.1.0D-6) go to 100
c
c
      write(6,*) 'iter=',iter,'time=',tsum
c
      if (iter.eq.10) go to 31
      go to 99
   31 call graphic
c
      stop
      end
c
      subroutine graphic
      implicit double precision (a-h,o-z)
      parameter(if=21,jf=21,kf=21)
      common T(if,jf,kf)
      common dx,dy,dz
      dimension x(if),y(jf),z(kf)
c
      x(1)=0.0
      do i=2,if
      x(i)=x(i-1)+dx
      enddo
      y(1)=0.0
      do j=2,jf
      y(j)=y(j-1)+dy
      enddo
      z(1)=0.0
      do k=2,kf
      z(k)=z(k-1)+dz
      enddo
c
      open(11,file='CUBICP.DAT',status='unknown', blank='null')
      write(11,111) if,jf,kf
      write(11,112) (((real(x(i)),i=1,if),j=1,jf),k=1,kf)
      write(11,112) (((real(y(j)),i=1,if),j=1,jf),k=1,kf)
      write(11,112) (((real(z(k)),i=1,if),j=1,jf),k=1,kf)
      close(11)
      open(11,file='DCUBICP.DAT',status='unknown', blank='null')
       uin=1.0
c       vis=0.0
c       tin=0.0
c       time=0.0
       write(11,111) if,jf,kf
       write(11,112) uin
c       write(11,112) (((1.0,i=1,if),j=1,jf),k=1,kf)
c       write(11,112) (((0.0,i=1,if),j=1,jf),k=1,kf)
c       write(11,112) (((0.0,i=1,if),j=1,jf),k=1,kf)
c       write(11,112) (((0.0,i=1,if),j=1,jf),k=1,kf)
       write(11,112) (((T(i,j,k),i=1,if),j=1,jf),k=1,kf)
      close(11)
c
  111 format(7I9)
  112 format(7F10.5)
c
      return
      end