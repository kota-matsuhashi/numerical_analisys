module mdl_schs_source
    use mdl_input_manager
    use mdl_mpisub_sbsp
    use mdl_decompo
    use mdl_param
    use mdl_refs
    use mdl_block
    use mdl_fdm1
    use mdl_halo
    use mdl_steps
    use mdl_table
    implicit none
  
    real(8), public, parameter :: SCHS1 = 1.0d3
    integer :: idata1=11
    integer :: idata2=7
  
    real(8), save, public, allocatable, dimension(:) :: SCHS_PKT1
    real(8), save, public, allocatable, dimension(:) :: SCHS_PKK1
    real(8), save, public, allocatable, dimension(:) :: SCHS_SLO1
    real(8), save, public, allocatable, dimension(:) :: SCHS_PKT2
    real(8), save, public, allocatable, dimension(:) :: SCHS_PKK2
    real(8), save, public, allocatable, dimension(:) :: SCHS_SLO2
    
  contains
  !-------------------------------------------------------------------------------
    subroutine schs_source(BLK,RHS)
      type(block), intent(inout) :: BLK
      real(8), intent(inout), dimension(:,:,:,:) :: RHS
      !real(8), allocatable, dimension(:,:,:) :: SOL
      
      ! integer :: i,j,k,imax,jmax,kmax
      ! integer :: nstep,step_start,step_final
      integer :: i,j,k
      integer :: ivalid
      integer :: myrank,id_block,nstep
      real(8) :: cfdim, srct, cklno, dckin
      integer :: n, iFLG_SCHS
      real(8) :: logkw,logkh,react
      real(8),parameter :: logrf=86.0d0,gasco=8.3144626d0,psi=-450000d0
      real(8) :: k_ref,k_int,ohmns,t_ref
      real(8) :: a1_ne,b1_ne,c1_ne,d1_ne,e1_ne,f1_ne,g1_ne
      real(8) :: a2_ne,b2_ne,c2_ne,d2_ne,e2_ne,f2_ne,g2_ne
      real(8) :: a3_ne,b3_ne,c3_ne,d3_ne
      real(8) :: kb,na,pi,sgma,gcon,rhni,v_nio,r_ave,C_solid,coag,e
      real(8) :: mo,ms,an,bt,rh,tp
      real(8) :: sups,c_mo,crdi,nuge,nugr,ni,izer,diff,logsups,solb,mass
  
      call mpisub_sbsp_get_myrank_world(myrank)
      call steps_get_step_now(nstep)
      
      100 format(100e14.6)
  
      call look_up_table_calc_func(TBL_TP_to_DIELEC,BLK%TMPR,TTREF,BLK%PRSS,PPREF, &
      BLK%DIEL,1.0d0,BLK%ista,BLK%iend,BLK%jsta,BLK%jend,BLK%ksta,BLK%kend,ivalid)
      
  
      a1_ne = -4.098d0
      b1_ne = -3245.2d0
      c1_ne = 223630.0d0
      d1_ne = -39984000.0d0
      e1_ne = 13.957d0
      f1_ne = -1262.3d0
      g1_ne = 856410.0d0
  
      a2_ne = -4.8d0
      b2_ne = 22.4d0
      c2_ne = -14.8d0
  
      a3_ne = -139.5d0
      b3_ne = -4508.74d0
      c3_ne = 14.26638d0
      d3_ne = 5.556747d0
  
      k_ref=86.0d0
      t_ref=380d0+273.15d0
      
      kb = 1.38065d-23   
      na = 6.02214d23  
      ni = 74.69d-3
      rhni = 6.67d3
      gcon = 8.3144626d0
      v_nio = ni/(rhni*na)
      pi = 4.0d0*DATAN(1.0d0)
      C_solid = rhni/ni
      e = DEXP(1.0d0)
       
      !! FURUKAWA
      !! MEMO--------------
      !! This roop is not vectorized due to the inner loop. Need the implementation for vectorizing.
      do k=BLK%ksta,BLK%kend
        do j=BLK%jsta,BLK%jend
          do i=BLK%ista,BLK%iend
            rh = BLK%DNST(i,j,k)*RHREF
            tp = BLK%TMPR(i,j,k)*TTREF
            ms = BLK%YSPC(i,j,k,1) 
            mo = BLK%YSPC(i,j,k,2)
            an = BLK%YSPC(i,j,k,3)
            bt = BLK%YSPC(i,j,k,4)
            
            solb = rh*e**(a3_ne-b3_ne/tp+c3_ne*DLOG(tp)+d3_ne*DLOG(rh))
                                              
            !mol/kgベースで計算
            logkh=a2_ne*(rh*1.0d-3)**2.0d0+b2_ne*rh*1.0d-3+c2_ne
            logkw=(a1_ne+b1_ne/tp+c1_ne/(tp**2.0d0)+d1_ne/(tp**3.0d0)) &
                  +(e1_ne+f1_ne/tp+g1_ne/(tp**2.0d0))*DLOG10(rh*1.0d-3)
            ohmns = 10**logkw/dsqrt(10**logkh)
            k_int = dexp((BLK%DIEL(i,j,k)-1.0d0)/(2.0d0*BLK%DIEL(i,j,k)+1)*psi/(gasco*t_ref)+k_ref)
            BLK%KNET(i,j,k)=k_int*(ohmns**2.0d0)
            srct  = BLK%KNET(i,j,k)*rh*ms
            cfdim = XYREF/(RHREF*UUREF)
            
  
            sups = rh*mo/solb !solb mol/m^3 
            diff = 2.24D-7*(tp)**0.763/(rh) 
            sgma = 0.1d0*kb*tp*(C_solid*na)**(2.0d0/3.0d0)*DLOG(C_solid/solb)
  
            crdi = 0.0d0
            nuge = 0.0d0
            nugr = 0.0d0
            if(sups>1.0d0) then
              crdi = 2.0d0*sgma*v_nio/(kb*tp*Dlog(sups))
              nuge = 1.5d0*diff*(na*mo*rh)**(7.0d0/3.0d0)*dsqrt(sgma/(kb*tp))*v_nio &
                    *dexp((-4.0d0*pi*crdi**2.0d0*sgma)/(3.0d0*kb*tp))
              nugr = 4.0d-5*(sups-1.0d0)
            end if
  
            nuge = DMIN1(1.0d30,nuge)
  
            BLK%SUPS(i,j,k) = sups
            BLK%NUGE(i,j,k) = nuge
            BLK%CRDI(i,j,k) = crdi 
  
            r_ave = 0.0d0
            if(bt>0.and.an>1.0d4) then
              r_ave = ((3.0d0*bt)/(4.0d0*pi*rhni*an))**(1.0d0/3.0d0)
            end if
  
            mass = (4.0d0*pi*crdi**3*nuge/3.0d0 + 4.0d0*pi*r_ave**2*nugr*rh*an)*rhni
            BLK%MASS(i,j,k) = mass
            BLK%SGMA(i,j,k) = r_ave
  
            coag = 3.5d-22
  
            !RHS(i,j,k,neqn_passive_scalar  ) = RHS(i,j,k,neqn_passive_scalar  ) - BLK%RJCB(i,j,k)*cfdim*srct
            !RHS(i,j,k,neqn_passive_scalar+1) = RHS(i,j,k,neqn_passive_scalar+1) + BLK%RJCB(i,j,k)*(srct &
                                                                -mass/ni)*cfdim
            !                                                             )*cfdim                                                                                                      
            !RHS(i,j,k,neqn_passive_scalar+2) = RHS(i,j,k,neqn_passive_scalar+2) + BLK%RJCB(i,j,k) &
                                  *(nuge-(rh*an)**2.0d0*coag)*cfdim
            !RHS(i,j,k,neqn_passive_scalar+3) = RHS(i,j,k,neqn_passive_scalar+3) + BLK%RJCB(i,j,k)*mass*cfdim
            RHS(i,j,k,neqn_passive_scalar  ) = RHS(i,j,k,neqn_passive_scalar  ) - BLK%RJCB(i,j,k)*cfdim*srct
            RHS(i,j,k,neqn_passive_scalar+1) = RHS(i,j,k,neqn_passive_scalar+1) + BLK%RJCB(i,j,k)*(cfdim*srct)
            RHS(i,j,k,neqn_passive_scalar+2) = RHS(i,j,k,neqn_passive_scalar+2)
            RHS(i,j,k,neqn_passive_scalar+3) = RHS(i,j,k,neqn_passive_scalar+3)
          enddo
        enddo
      enddo
  
      if(nstep==80000) then
        if(myrank==75) then
          open(300,file='data/writea.txt',form='formatted')
          write(300,*) "nstep=0"
          write(300,100) BLK%SUPS(BLK%ista+100,BLK%jsta+2,BLK%ksta+10)
          write(300,100) BLK%NUGE(BLK%ista+100,BLK%jsta+2,BLK%ksta+10)
          write(300,100) BLK%SGMA(BLK%ista+100,BLK%jsta+2,BLK%ksta+10)
          write(300,100) BLK%CRDI(BLK%ista+100,BLK%jsta+2,BLK%ksta+10)
          write(300,100) BLK%DNST(BLK%ista+100,BLK%jsta+2,BLK%ksta+10) *RHREF
          write(300,100) BLK%TMPR(BLK%ista+100,BLK%jsta+2,BLK%ksta+10) *TTREF
          write(300,100) BLK%YSPC(BLK%ista+100,BLK%jsta+2,BLK%ksta+10,2)
          close(300)
        end if
      end if
      if(nstep==80001) then
        if(myrank==75) then
          open(300,file='data/writeb.txt',form='formatted')
          write(300,*) "nstep=1"
          write(300,100) BLK%SUPS(BLK%ista+100,BLK%jsta+2,BLK%ksta+10)
          write(300,100) BLK%NUGE(BLK%ista+100,BLK%jsta+2,BLK%ksta+10)
          write(300,100) BLK%SGMA(BLK%ista+100,BLK%jsta+2,BLK%ksta+10)
          write(300,100) BLK%CRDI(BLK%ista+100,BLK%jsta+2,BLK%ksta+10)
          write(300,100) BLK%DNST(BLK%ista+100,BLK%jsta+2,BLK%ksta+10) *RHREF
          write(300,100) BLK%TMPR(BLK%ista+100,BLK%jsta+2,BLK%ksta+10) *TTREF
          write(300,100) BLK%YSPC(BLK%ista+100,BLK%jsta+2,BLK%ksta+10,2)
          close(300)
        end if
      end if
  
    end subroutine schs_source
  
    subroutine schs_init
      integer :: myrank,nprocs,iproc
      integer :: ierror,i
  #ifdef _USE_MPI
      call MPI_Comm_rank(MPI_COMM_WORLD,myrank,ierror)
      call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierror)
  #else
      myrank = 0
      nprocs = 1
  #endif
      allocate(SCHS_PKT1(idata1))
      allocate(SCHS_PKK1(idata1))
      allocate(SCHS_SLO1(idata1))
      allocate(SCHS_PKT2(idata2))
      allocate(SCHS_PKK2(idata2))
      allocate(SCHS_SLO2(idata2))
      
      do iproc=0,nprocs-1
        if (myrank==iproc) then
          open(11,file='./data/solubility/Ni_solubility.txt',form='formatted')
           read(11,*) (SCHS_PKT1(i),i=1,idata1)
           read(11,*) (SCHS_PKK1(i),i=1,idata1)
           read(11,*) (SCHS_SLO1(i),i=1,idata1)
        endif
      end do
  
      do iproc=0,nprocs-1
        if (myrank==iproc) then
          open(12,file='./data/reaction_rate/Ni.txt',form='formatted')
           read(12,*) (SCHS_PKT2(i),i=1,idata2)
           read(12,*) (SCHS_PKK2(i),i=1,idata2)
        endif
      end do
      
    end subroutine schs_init
    subroutine schs_final
      deallocate(SCHS_PKT1)
      deallocate(SCHS_PKK1)
      deallocate(SCHS_SLO1)
      deallocate(SCHS_PKT2)
      deallocate(SCHS_PKK2)
      deallocate(SCHS_SLO1)
    end subroutine schs_final
  
  end module mdl_schs_source