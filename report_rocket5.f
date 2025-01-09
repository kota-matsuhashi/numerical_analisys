C----------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter(if=121,jf=21,kf=45)
      common T(if,jf,kf)
      common dx,dy,dz
      dimension al(if,jf,kf)
      real :: time_begin_s,time_end_s
C
      data omega/1.5/
c---層の分離
      data ia,ib,ic/1,41,41/
      data ja,jb,jc/1,21,21/
      data ka,kb,kc/37,39,45/
c---熱伝導係数
c---1層目(空気)
      data ck1/26.14/
c---2層目(アブレータ)
      data ck2/0.42/
c---3層目(モータケース)
      data ck3/120/
c---初期値の設定
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
c---計算時間の計測
      call cpu_time(time_begin_s)
C---境界条件の設定
      do k=1,kf
      do j=1,jf
      do i=1,if
        T(i,j,k)=30.0
        al(i,j,k)=30.0
      enddo
      enddo
      enddo
C---i方向熱勾配(流路方向)
c---空気
      do k=1,ka - 1
      do j=1,jf
        T(1,j,k)=30.0
        T(if,j,k)=30.0
      enddo
      enddo
c---アブレータ,モータケース断面
      do k=ka,kf
      do j=1,jf
c        T(1,j,k)=500.0
        T(if,j,k)=30.0
      enddo
      enddo 
c---k方向熱勾配(半径方向)
      do i=1,if
      do j=1,jf
        T(i,j,1)=30.0
        T(i,j,kf)=30.0
      enddo
      enddo
c---イグナイタ
      T(1, 11, 1)=800.0
      T(1, 12, 1)=800.0
      T(1, 10, 1)=800.0
      T(2, 11, 1)=800.0
      T(2, 12, 1)=800.0
      T(2, 10, 1)=800.0 
      T(1, 11, 2)=800.0
      T(1, 12, 2)=800.0
      T(1, 10, 2)=800.0   
C---計算時間インクリメント
      iter=0
      tsum=0.0
   99 iter=iter+1
      tsum=tsum+dt
c---alfaの計算---------------------------
      do i=2,if-1
      do j=1,jf
      do k=1,kf
c---x方向の熱伝導係数
      im1=i-1
C     if (i.eq.1) im1=2 doループで除外
      ip1=i+1
C     if (i.eq.if) ip1=if-1 doループで除外
c---1層目
      ckxl=ck1
      ckxr=ck1
c---2層目
      if ((j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (i.gt.ia).and.(i.le.ib)) ckxl=ck2
      if ((j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (i.ge.ia).and.(i.lt.ib)) ckxr=ck2
c---3層目
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (j.gt.jb).and.(j.le.jc)) ckxl=ck3
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (j.ge.jb).and.(j.lt.jc)) ckxr=ck3
c---y方向の熱伝導係数
      jm1=j-1
      if (j.eq.1) jm1=2
      jp1=j+1
      if (j.eq.jf) jp1=jf-1
c---1層目
      ckyl=ck1
      ckyr=ck1
c---2層目
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (j.gt.ja).and.(j.le.jb)) ckyl=ck2
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (j.ge.ja).and.(j.lt.jb)) ckyr=ck2
c---3層目
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (j.gt.jb).and.(j.le.jc)) ckyl=ck3
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (j.ge.jb).and.(j.lt.jc)) ckyr=ck3
c---z方向の熱伝導係数
      km1=k-1
      if (k.eq.1) km1=2
      kp1=k+1
      if (k.eq.kf) kp1=kf-1
c---1層目
      ckzl=ck1
      ckzr=ck1
c---2層目
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (j.ge.ja).and.(j.le.jb).and.
     &    (k.gt.ka).and.(k.le.kb)) ckzl=ck2
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.lt.kb)) ckzr=ck2
c---3層目
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (j.gt.jb).and.(j.le.jc)) ckzl=ck3
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (j.ge.jb).and.(j.lt.jc)) ckzr=ck3
C---alpa
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
C---メイン計算ループ-----------------------------------
      na=0
  100 na=na+1
      rmax=0.0
c
      do ii=2,if-1
      i=ii
      if (mod(na,2).eq.0) i=if+1-i
      do jj=1,jf
      j=jj
      if (mod(na,2).eq.0) j=jf+1-j
      do kk=1,kf
      k=kk
      if (mod(na,2).eq.0) k=kf+1-k
c---x方向の熱伝導係数
      im1=i-1
C     if (i.eq.1) im1=2
      ip1=i+1
C     if (i.eq.if) ip1=if-1
c---1層目
      ckxl=ck1
      ckxr=ck1
c---2層目
      if ((j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (i.gt.ia).and.(i.le.ib)) ckxl=ck2
      if ((j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (i.ge.ia).and.(i.lt.ib)) ckxr=ck2
c---3層目
      if ((j.ge.jb).and.(j.le.jc).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (i.gt.ib).and.(i.le.ic)) ckxl=ck3
      if ((j.ge.jb).and.(j.le.jc).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (i.ge.ib).and.(i.lt.ic)) ckxr=ck3
c---y方向の熱伝導係数
      jm1=j-1
      if (j.eq.1) jm1=2
      jp1=j+1
      if (j.eq.jf) jp1=jf-1
c---1層目
      ckyl=ck1
      ckyr=ck1
c---2層目
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (j.gt.ja).and.(j.le.jb)) ckyl=ck2
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (k.ge.ka).and.(k.le.kb).and.
     &    (j.ge.ja).and.(j.lt.jb)) ckyr=ck2
c---3層目
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (j.gt.jb).and.(j.le.jc)) ckyl=ck3
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (k.ge.kb).and.(k.le.kc).and.
     &    (j.ge.jb).and.(j.lt.jc)) ckyr=ck3
c---z方向の熱伝導係数
      km1=k-1
      if (k.eq.1) km1=2
      kp1=k+1
      if (k.eq.kf) kp1=kf-1
c---1層目
      ckzl=ck1
      ckzr=ck1
c---2層目
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (j.ge.ja).and.(j.le.jb).and.
     &    (k.gt.ka).and.(k.le.kb)) ckzl=ck2
      if ((i.ge.ia).and.(i.le.ib).and.
     &    (j.ge.ja).and.(j.le.jb).and.
     &    (k.ge.ka).and.(k.lt.kb)) ckzr=ck2
c---3層目
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (j.ge.jb).and.(j.le.jc).and.
     &    (k.gt.kb).and.(k.le.kc)) ckzl=ck3
      if ((i.ge.ib).and.(i.le.ic).and.
     &    (j.ge.jb).and.(j.le.jc).and.
     &    (k.ge.kb).and.(k.lt.kc)) ckzr=ck3
c---Laplace演算子の計算
      dl=1.0+(dlx*(ckxr+ckxl)
     &       +dly*(ckyr+ckyl)
     &       +dlz*(ckzr+ckzl))*0.5
c---r=T(m+1)-T(m)
      r=(ckxl*T(im1,j,k)+ckxr*T(ip1,j,k))*dlx/(dl*2.0)
     & +(ckyl*T(i,jm1,k)+ckyr*T(i,jp1,k))*dly/(dl*2.0)
     & +(ckzl*T(i,j,km1)+ckzr*T(i,j,kp1))*dlz/(dl*2.0)
     &  +al(i,j,k)/dl - T(i,j,k)
c---T(m+1)
      T(i,j,k)=T(i,j,k) + omega*r
c---rmaxを正の数に変換
      rmax=amax1(rmax,abs(r))
      enddo
      enddo
      enddo
c---rmaxを正規化
      if (na.eq.1) rmax1=rmax
      rmax=rmax/rmax1
      write(6,20) na, rmax
      write(11,20) na,rmax
   20 format(5x,'na=',i3,5x,'rmax=',f10.7)
c---収束判定
      if (na.lt.100.and.rmax.gt.1.0D-6) go to 100
c---収束した場合
      write(6,*) 'iter=',iter,'time=',tsum
c---iteration指定
      if (iter.eq.1000) go to 31
c---iter=10未満の場合
      go to 99
   31 call graphic
c---計算時間の計測
      call cpu_time(time_end_s)
      print *,time_begin_s, time_end_s
      print *,time_end_s - time_begin_s,"sec"
      stop
      end
c---データ出力サブルーチン
      subroutine graphic
      implicit double precision (a-h,o-z)
      parameter(if=121,jf=21,kf=45)
      common T(if,jf,kf)
      common dx,dy,dz
      dimension x(if),y(jf),z(kf)
c---グリッドデータの作成
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
c---グリッドデータの書き出し
      open(11,file='grid_rocket.DAT',status='unknown', blank='null')
      write(11,111) if,jf,kf
      write(11,112) (((real(x(i)),i=1,if),j=1,jf),k=1,kf)
      write(11,112) (((real(y(j)),i=1,if),j=1,jf),k=1,kf)
      write(11,112) (((real(z(k)),i=1,if),j=1,jf),k=1,kf)
      close(11)
c---データファイルの書き出し
      open(11,file='result_rocket5.DAT',status='unknown', blank='null')
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
       write(11,112) ((( T(i,j,k),i=1,if),j=1,jf),k=1,kf)
      close(11)
c
  111 format(7I9)
  112 format(7F10.5)
c
      return
      end
