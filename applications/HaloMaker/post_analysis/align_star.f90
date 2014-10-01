!!$ from S. Courty
!!$ revised version by Leo (feb 2012)
!!$ - run from rats script (arguments added)
!!$ - to fit in memory limits, align the small sub-volume (subvol_...)
program align_star

  use ramses_util
  implicit none

  real(kind=8), parameter :: Msol = 1.9891d+33
  real(kind=8), parameter :: yr = 3.15569d+07
  real(kind=8), parameter :: kparsec = 3.08d+21
  real(kind=8), parameter :: kms = 1.d+05

  integer :: i,ii,j,ic,choice,n,nb,nstep,nfile,st
  integer::nstar
  integer::n_frw
  integer,dimension(:),allocatable::idstar,id,idj
  real(kind=8) :: aexp,boxlen,scale_l,scale_d,scale_t,scale_v
  real(kind=8) :: scale_cosmict,xc,omega_m,omega_l,omega_k,t,tau
  real(kind=8) :: amass,aexp_ini,time_tot,time_simu,ageu
  real(kind=8),dimension(:,:),allocatable :: xstar,vstar,xgas,vgas,pos,vel
  real(kind=8),dimension(:),allocatable :: astar,zstar,mstar,astar2,mass,met,ap,ap2
  real(kind=8),dimension(:),allocatable :: mgas,zgas,rhop,pp,dl,tp
  real(kind=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  character(len=256) :: arg
  character(len=256) :: repository,foutstar,fcell,fname,fstar
#ifdef CHEMO
  integer :: nelt
  real(kind=8),dimension(:),allocatable :: tstar,m0star,mass0,ap3
  real(kind=8),dimension(:,:),allocatable :: eltstar,elt
#endif

  integer :: nleaf,lmax,nvarh,iimax
  real(kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax      
  real(kind=8) :: Js,Jxs,Jys,Jzs,mtot,jjx,jjy,jjz,limit,size, &
       ax,bx,cx,ay,by,cy,az,bz,cz,txs,tys,tzs, &
       dx,dy,dd,costh,sinth,cosph,sinph,d1,d2,r0,xc2

  real(kind=4) :: rlim,rji
  integer :: iargc
  character(len=4)   :: opt
  character(len=256) :: outdir,subdir
  logical :: verbose

  !!======================================================================================
  !! Method to compute angular momentum
  !! 0 --> gas
  !! 1 --> stars
  choice=0

  !! read params
  !! set default values
  rji  = 5.
  outdir = "."
  subdir = "."
  verbose = .false.
  n = iargc()
  if (n < 2) then
     print *, 'usage: align_star  -nst nstep'
     print *, '                   [-sub subvol_dir] '
     print *, '                   [-out output_dir]'
     print *, '                   [-rji rlim_angmom] '
     print *, '                   [-ver verbose] '
     stop
  end if
  do i=1,n,2
     call getarg(i,opt)
     if (i == n) then
        print '("option ",a4," has no argument")', opt
        stop
     end if
     call getarg(i+1,arg)
     select case (opt)
     case ('-nst')
        read(arg,*) nfile
     case ('-out')
        outdir = trim(arg)
     case ('-sub')
        subdir = trim(arg)
     case ('-rji')
        read(arg,*) rji
     case ('-ver')
        verbose=.true.
     case default
        print '("unknown option ",a4," ignored")', opt
     end select
  end do

  !!======================================================================================
  print *
  write(arg,'("/subvol_star_00",I3,".dat")') nfile
  if (nfile<100) write(arg,'("/subvol_star_000",I2,".dat")') nfile
  if (nfile<10) write(arg,'("/subvol_star_0000",I1,".dat")') nfile
  fstar = trim(subdir)//trim(arg)
  print *,'read file: ',trim(fstar)
  open(13, file=fstar, status='unknown', form='unformatted')
  rewind 13
  read(13) r0,scale_l,scale_d,scale_t,aexp,ageu,time_simu
  read(13) nstar
  print *,'nstar =',nstar
  allocate(xstar(nstar,3),vstar(nstar,3),astar(nstar),mstar(nstar),astar2(nstar),idstar(nstar),zstar(nstar))
  read(13) xstar(:nstar,:3)
  read(13) vstar(:nstar,:3)
  read(13) mstar(:nstar)
  read(13) astar(:nstar)
  read(13) astar2(:nstar)
  read(13) zstar(:nstar)
  read(13) idstar(:nstar)
#ifdef CHEMO
     allocate(m0star(nstar),tstar(nstar))
     read(13) m0star(:nstar)
     read(13) tstar(:nstar)
     read(13) nelt
     allocate(eltstar(nstar,nelt))
     read(13) eltstar(:nstar,:nelt)
#endif
  close(13)

  if (choice==0) then  

     write(arg,'("/subvol_gas_00",I3,".dat")') nfile
     if (nfile<100) write(arg,'("/subvol_gas_000",I2,".dat")') nfile
     if (nfile<10) write(arg,'("/subvol_gas_0000",I1,".dat")') nfile
     fcell = trim(subdir)//trim(arg)
     print *,'read file: ',trim(fcell)
     open(unit=10,file=fcell,form='unformatted')
     read(10) r0,scale_l,scale_d,scale_t,aexp
     read(10) nleaf
     print *,'nleaf =',nleaf
     allocate(xgas(nleaf,3),vgas(nleaf,3),mgas(nleaf),zgas(nleaf),rhop(nleaf),pp(nleaf),dl(nleaf),tp(nleaf))
     read(10) xgas(:nleaf,:3)
     read(10) vgas(:nleaf,:3)
     read(10) mgas(:nleaf)
     read(10) zgas(:nleaf)
     read(10) rhop(:nleaf)
     read(10) pp(:nleaf)
     read(10) dl(:nleaf)
     read(10) tp(:nleaf)
     close(10) 

  endif

  !!======================================================================================
  amass=scale_d*scale_l**3
  amass=amass/Msol
  print '(a8,1pe12.5)','amass =',amass

  size=scale_l/3.08d21
  print *,'size (kpc) =',size
  print *

  !!======================================================================================
  !! Leo: now align stars

  print *,'align data according to the angular momentum computed in a sphere of radius'

  limit=rji/size
  print '(a7,2x,e11.4,2x,f6.2)','rlim =',limit,rji

  ii=0
  mtot=0.
  Jxs=0.
  Jys=0.
  Jzs=0.
  Js=0.

  if (choice==0) then  
     print *,'...using gas'
     do j=1,nleaf 
        dd=sqrt(xgas(j,1)*xgas(j,1) + xgas(j,2)*xgas(j,2) + xgas(j,3)*xgas(j,3))
        if (dd<limit) then
           ii=ii+1
           mtot = mtot+mgas(j)
           Jxs=Jxs+mgas(j)*(xgas(j,2)*vgas(j,3) - xgas(j,3)*vgas(j,2))
           Jys=Jys+mgas(j)*(xgas(j,3)*vgas(j,1) - xgas(j,1)*vgas(j,3))
           Jzs=Jzs+mgas(j)*(xgas(j,1)*vgas(j,2) - xgas(j,2)*vgas(j,1))            
        endif
     enddo
     deallocate(xgas,vgas,mgas,zgas,rhop,pp,dl,tp)
  else 
     print *,'...using stars'
     do j=1,nstar
        dd=sqrt(xstar(j,1)*xstar(j,1) + xstar(j,2)*xstar(j,2) + xstar(j,3)*xstar(j,3))
        !???????
        if ( (dd<limit) .and. (astar2(j)<12.5)) then
           ii=ii+1
           mtot = mtot+mstar(j)
           Jxs=Jxs+mstar(j)*(xstar(j,2)*vstar(j,3) - xstar(j,3)*vstar(j,2))
           Jys=Jys+mstar(j)*(xstar(j,3)*vstar(j,1) - xstar(j,1)*vstar(j,3))
           Jzs=Jzs+mstar(j)*(xstar(j,1)*vstar(j,2) - xstar(j,2)*vstar(j,1))
        endif
     enddo
  endif
  iimax=ii
  print *,'iimax =',iimax
  print *,'mtot =',mtot
  Js=sqrt(Jxs*Jxs+Jys*Jys+Jzs*Jzs)
  print *,'Js ',Js,Jxs,Jys,Jzs


  if (Js > 0.) then 
     jjx=Jxs/Js
     jjy=Jys/Js
     jjz=Jzs/Js
     costh=jjz
     sinth=sqrt(1.0-jjz*jjz)
     if  (sinth > 0.0) then 
        sinph=jjy/sinth
        cosph=jjx/sinth
     endif
  endif
  if (Js <= 0.)  then 
     cosph = 1.0
     sinph = 0.0
  endif

  ax=costh*cosph
  bx=costh*sinph
  cx=-sinth
  ay=-sinph
  by=cosph
  cy=0.0
  az=sinth*cosph
  bz=sinth*sinph
  cz=costh

  print *,'ax,bx,cx'
  print *, ax,bx,cx
  print *, ay,by,cy
  print *, az,bz,cz

  do j=1,nstar
     txs=xstar(j,1)
     tys=xstar(j,2)
     tzs=xstar(j,3)
     xstar(j,1)=(ax*txs+bx*tys+cx*tzs)
     xstar(j,2)=(ay*txs+by*tys+cy*tzs)
     xstar(j,3)=(az*txs+bz*tys+cz*tzs)

     txs=vstar(j,1)
     tys=vstar(j,2)
     tzs=vstar(j,3)
     vstar(j,1)=(ax*txs+bx*tys+cx*tzs)
     vstar(j,2)=(ay*txs+by*tys+cy*tzs)
     vstar(j,3)=(az*txs+bz*tys+cz*tzs)         
  enddo

  !!  scale_v=scale_l/scale_t/kms
  !!  vstar(:nstar,:3)=vstar(:nstar,:3)*scale_v

  !!======================================================================================
  !! Leo: select a sphere and dump align file

  ! r0 read in subvol header
  allocate(idj(:nstar))
  ii=0
  idj=0
  do j=1,nstar
     dd=xstar(j,1)*xstar(j,1)+xstar(j,2)*xstar(j,2)+xstar(j,3)*xstar(j,3)
     if (sqrt(dd)>r0) cycle
     ii=ii+1           
     idj(ii) = j
  enddo
  iimax=ii
  print *,'iimax =',iimax

  allocate(pos(:iimax,:3),vel(:iimax,:3),mass(:iimax),ap(:iimax),ap2(:iimax),met(:iimax),id(:iimax))
#ifdef CHEMO
  allocate(elt(:iimax,:nelt),mass0(:iimax),ap3(:iimax))
#endif
  do ii=1,iimax
     j = idj(ii)
     pos(ii,1)=xstar(j,1)
     pos(ii,2)=xstar(j,2)
     pos(ii,3)=xstar(j,3)
     vel(ii,1)=vstar(j,1)
     vel(ii,2)=vstar(j,2)
     vel(ii,3)=vstar(j,3)
     mass(ii)=mstar(j)
     ap(ii)=astar(j)
     ap2(ii)=astar2(j)
     met(ii)=zstar(j)
     id(ii)=idstar(j)
#ifdef CHEMO     
     mass0(ii)=m0star(j)
     ap3(ii)=tstar(j)
     elt(ii,:nelt)=eltstar(j,:nelt)
#endif
  enddo
  deallocate(xstar,mstar,vstar,idstar,zstar,astar,astar2)
  deallocate(idj)
#ifdef CHEMO
  deallocate(m0star,tstar,eltstar)
#endif
  nstar=iimax

  write(arg,'("/align_star_00",I3,".dat")') nfile
  if (nfile<100) write(arg,'("/align_star_000",I2,".dat")') nfile
  if (nfile<10) write(arg,'("/align_star_0000",I1,".dat")') nfile
  foutstar = trim(outdir)//trim(arg)

  open(13, file=foutstar, status='unknown', form='unformatted') 
  rewind 13 
  write(13) r0,scale_l,scale_d,scale_t,aexp,ageu,time_simu
  write(13) nstar
  write(13) pos(:nstar,:3)
  write(13) vel(:nstar,:3)
  write(13) mass(:nstar)
  write(13) ap(:nstar)
  write(13) ap2(:nstar)
  write(13) met(:nstar)
  write(13) id(:nstar)
#ifdef CHEMO
  write(13) mass0(:nstar)
  write(13) ap3(:nstar)
  write(13) nelt
  write(13) elt(:nstar,:nelt)
#endif
  close(13)

  print *,'=> dump file: ',trim(foutstar)
  print *
  deallocate(pos,vel,mass,ap,ap2,met,id)
#ifdef CHEMO
  deallocate(ap3,mass0,elt)
#endif

end program align_star
