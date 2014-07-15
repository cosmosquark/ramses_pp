!!$ from S. Courty
!!$ revised version by Leo (feb 2012)
!!$ - run from rats script (arguments added)
!!$ - to fit in memory limits, align the small sub-volume (subvol_...)
program align_dm

  use ramses_util
  implicit none

  real(kind=8), parameter :: Msol = 1.9891d+33
  real(kind=8), parameter :: yr = 3.15569d+07
  real(kind=8), parameter :: kparsec = 3.08d+21
  real(kind=8), parameter :: kms = 1.d+05

  integer :: i,ii,j,ic,choice,n,nb,nstep,nfile,st
  integer :: nstar,npart
  integer :: n_frw
  integer,dimension(:),allocatable::id,idstar,idp,idj
  real(kind=8) :: aexp,boxlen,scale_l,scale_d,scale_t,scale_v,mpdm
  real(kind=8) :: scale_cosmict,zz,xc,omega_m,omega_l,omega_k,t,tau
  real(kind=8) :: scale,amass,aexp_ini,time_tot,time_simu,ageu
  real(kind=8),dimension(:,:),allocatable::xp,vp,xg,xstar,vstar,xgas,vgas,pos,vel
  real(kind=8),dimension(:),allocatable:: astar,zstar,astar2,mstar,mgas,zgas,rhop,dl,pp,tp
  real(kind=8),dimension(:),allocatable:: aexp_frw,hexp_frw,tau_frw,t_frw
  character(len=256) :: arg
  character(len=256) :: repository,fcell,foutdm,fname,fdm,fstar

  integer :: nleaf,lmax,nvarh,iimax
  real(kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax      
  real(kind=8) :: Js,Jxs,Jys,Jzs,mtot,jjx,jjy,jjz,limit,size, &
       ax,bx,cx,ay,by,cy,az,bz,cz,txs,tys,tzs, &
       dx,dy,dd,costh,sinth,cosph,sinph,d1,d2,r0,xc2

  real(kind=4) :: rji
  integer       :: iargc
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
     print *, 'usage: align_dm   -nst nstep'
     print *, '                 [-sub subvol_dir] '
     print *, '                 [-out output_dir]'
     print *, '                 [-rji rlim_angmom] '
     print *, '                 [-ver verbose] '
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
     case ('-out')
        outdir = trim(arg)
     case ('-sub')
        subdir = trim(arg)
     case ('-nst')
        read(arg,*) nfile
     case ('-rji')
        read(arg,*) rji
     case ('-ver')
        verbose=.true.
     case default
        print '("unknown option ",a4," ignored")', opt
     end select
  end do

  !!======================================================================================

  write(arg,'("/subvol_dm_00",I3,".dat")') nfile
  if (nfile<100) write(arg,'("/subvol_dm_000",I2,".dat")') nfile
  if (nfile<10) write(arg,'("/subvol_dm_0000",I1,".dat")') nfile
  fdm = trim(subdir)//trim(arg)
  print *,'read file: ',trim(fdm)
  open(13, file=fdm, status='unknown', form='unformatted')
  rewind 13

  read(13) r0,scale_l,scale_d,scale_t,aexp,mpdm
  read(13) npart
  print *,'npart =',npart
  allocate(xp(npart,3),vp(npart,3),idp(npart))
  read(13) xp(:npart,:3)
  read(13) vp(:npart,:3)
  read(13) idp(:npart)
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

  else

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
     close(13)

  endif

  print *

  !!======================================================================================

  amass=scale_d*scale_l**3
  amass=amass/Msol
  print '(a8,1pe12.5)','amass =',amass

  size=scale_l/3.08d21
  print *,'size (kpc) =',size
  print *

  !!======================================================================================

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
     deallocate(xstar,vstar,mstar,zstar,astar,astar2)
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

  do j=1,npart
     txs=xp(j,1)
     tys=xp(j,2)
     tzs=xp(j,3)
     xp(j,1)=(ax*txs+bx*tys+cx*tzs)
     xp(j,2)=(ay*txs+by*tys+cy*tzs)
     xp(j,3)=(az*txs+bz*tys+cz*tzs)

     txs=vp(j,1)
     tys=vp(j,2)
     tzs=vp(j,3)
     vp(j,1)=(ax*txs+bx*tys+cx*tzs)
     vp(j,2)=(ay*txs+by*tys+cy*tzs)
     vp(j,3)=(az*txs+bz*tys+cz*tzs)         
  enddo

!!  scale_v=scale_l/scale_t/kms
!!  vp(:npart,:3)=vp(:npart,:3)*scale_v

  !!======================================================================================
  !! Leo: select a sphere and dump align file

  ! r0 read in subvol header
  allocate(idj(:npart))
  ii=0
  idj=0
  do j=1,npart
     dd=xp(j,1)*xp(j,1)+xp(j,2)*xp(j,2)+xp(j,3)*xp(j,3)
     if (sqrt(dd)>r0) cycle
     ii=ii+1           
     idj(ii) = j
  enddo
  iimax=ii
  print *,'iimax =',iimax

  allocate(pos(:iimax,:3),vel(:iimax,:3),id(:iimax))
  do ii=1,iimax
     j = idj(ii)
     pos(ii,1)=xp(j,1)
     pos(ii,2)=xp(j,2)
     pos(ii,3)=xp(j,3)
     vel(ii,1)=vp(j,1)
     vel(ii,2)=vp(j,2)
     vel(ii,3)=vp(j,3)
     id(ii)=idp(j)
  enddo
  deallocate(xp,vp,idp)
  deallocate(idj)

  npart=iimax

  write(arg,'("/align_dm_00",I3,".dat")') nfile
  if (nfile<100) write(arg,'("/align_dm_000",I2,".dat")') nfile
  if (nfile<10) write(arg,'("/align_dm_0000",I1,".dat")') nfile
  foutdm = trim(outdir)//trim(arg)

  open(13, file=foutdm, status='unknown', form='unformatted') 
  rewind 13 
  write(13) r0,scale_l,scale_d,scale_t,aexp,mpdm
  write(13) npart
  write(13) pos(:npart,:3)
  write(13) vel(:npart,:3)
  write(13) id(:npart)
  close(13) 
  print *,'=> dump file: ',trim(foutdm)
  print *
  deallocate(pos,vel,id)

end program align_dm
