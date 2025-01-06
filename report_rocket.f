C----------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(nr=21,ntheta=21,nz=21)
      common T(nr,ntheta,nz)
      !common /grid_params/ dr,dtheta,dz
      real :: dr, dtheta, dz
      dimension al(nr,ntheta,nz)
C
      data omega/1.5/
      data pi/3.141592653589793/
      data ia,ib/8,14/
      data ja,jb/8,14/
      data ka,kb/8,14/
      data ck1,ck2/10.0,1.0/
C-------- 格子間隔設定 --------
      dt=0.001
      r_inner = 0.5
      r_outer = 2.0
      dr=(r_outer - r_inner) / (nr - 1)
      dtheta=2.0*pi/ntheta
      dz=0.1
      dr2=dr*dr
      dtheta2=dtheta*dtheta
      dz2=dz*dz
      dlr=dt/dr2
      dltheta=dt/dtheta2
      dlz=dt/dz2
C------初期条件設定---------------------------------------------
      do k=1,nz
      do j=1,ntheta
      do i=1,nr
        T(i,j,k)=50.0
        al(i,j,k)=50.0
      enddo
      enddo
      enddo
C-----温度勾配を設定--------------------------------------------
      do k=1,nz
      do j=1,ntheta
        T(1,j,k)=100.0
        T(nr,j,k)=50.0
      enddo
      enddo
C
      iter=0
      tsum=0.0
   99 iter=iter+1
      tsum=tsum+dt
      do i=2,nr-1
      do j=1,ntheta
      do k=1,nz
c
      im1=i-1
C     if (i.eq.1) im1=2
      ip1=i+1
C     if (i.eq.nr) ip1=nr-1
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
      if (j.eq.ntheta) jp1=ntheta-1
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
      if (k.eq.nz) kp1=nz-1
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
c--------反復計算ルーチン--------------------------------------
      na=0
  100 na=na+1
      rmax=0.0
C
      do ii=2,nr-1
      i=ii
      if (mod(na,2).eq.0) i=nr+1-i
      do jj=1,ntheta
      j=jj
      if (mod(na,2).eq.0) j=ntheta+1-j
      do kk=1,nz
      k=kk
      if (mod(na,2).eq.0) k=nz+1-k
c
      im1=i-1
C     if (i.eq.1) im1=2
      ip1=i+1
C     if (i.eq.nr) ip1=nr-1
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
      if (j.eq.ntheta) jp1=ntheta-1
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
      if (k.eq.nz) kp1=nz-1
      ckzl=ck1
      ckzr=ck1
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (j.ge.ja).and.(j.le.jb).and.
     &    (k.gt.ka).and.(k.le.kb)) ckzl=ck2
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.lt.kb)) ckzr=ck2
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
      print *, 'dr=', dr
      print *, 'dtheta=', dtheta
      print *, 'dz=', dz
c
      if (iter.eq.10) go to 31
      go to 99
   31 call graphic
c
      stop
      end
c---------グラフィック表示サブルーチン---------------------------
      subroutine graphic
      implicit double precision (a-h,o-z)
      parameter(nr=21,ntheta=21,nz=21)
      common T(nr,ntheta,nz), al(nr,ntheta,nz)
      !common /grid_params/ dr,dtheta,dz
      dimension r(nr),theta(ntheta),z(nz)
      integer :: i, j, k
      real :: dr, dtheta, dz
      data pi/3.141592653589793/

C-------- 格子点の座標を生成 --------
      dr=0.1
      r_inner = 0.5
      r_outer = 2.0
      dr=(r_outer - r_inner) / (nr - 1)
      dtheta=2.0*pi/ntheta
      dz=0.1
      print *, 'dr=', dr
      print *, 'dtheta=', dtheta
      print *, 'dz=', dz
      r(1) = r_inner
      do i = 2, nr
        r(i) = r(i-1) + dr
      end do
      theta(1) = 0.0
      do j = 2, ntheta
        theta(j) = theta(j-1) + dtheta
      end do
      z(1) = 0.0
      do k = 2, nz
        z(k) = z(k-1) + dz
      end do

C-------- 格子データの出力 --------
      open(unit=10, file='cylindrical_grid_data.DAT', status='unknown')
      write(10,*) nz, ntheta, nr
      do k = 1, nz
        do j = 1, ntheta
          do i = 1, nr
            write(10,*) r(i), theta(j), z(k)
            !write(10,*) i, j, k, r(i), theta(j), z(k), T(i, j, k)
            !&, al(i, j, k)
          end do
        end do
      end do
      close(10)
      open(11,file='grid_data2.DAT',status='unknown', blank='null')
      write(11,111) nr, ntheta, nz
      write(11,112) (((real(r(i)),i=1,nr),j=1,ntheta),k=1,nz)
      write(11,112) (((real(theta(j)),i=1,nr),j=1,ntheta),k=1,nz)
      write(11,112) (((real(z(k)),i=1,nr),j=1,ntheta),k=1,nz)
      close(11)
c------データファイルへの書き出し-------------------------------
      open(11,file='result_rocket.DAT',status='unknown', blank='null')
       uin=1.0
       vis=0.0
       tin=0.0
       time=0.0
       write(11,111) nr,ntheta,nz
       write(11,112) uin,vis,tin,time
       !write(11,112) (((1.0,i=1,nr),j=1,ntheta),k=1,nz)
       !write(11,112) (((0.0,i=1,nr),j=1,ntheta),k=1,nz)
       !write(11,112) (((0.0,i=1,nr),j=1,ntheta),k=1,nz)
       !write(11,112) (((0.0,i=1,nr),j=1,ntheta),k=1,nz)
       write(11,112) (((T(i,j,k),i=1,nr),j=1,ntheta),k=1,nz)
      close(11)
c
  111 format(7I9)
  112 format(7F10.5)
c
      return
      end