!!$ from S. Courty
!!$ revised version by Leo (feb 2012)
!!$ - run from rats script (arguments added)
!!$ - to fit in memory limits, align the small sub-volume (subvol_...)
program align_gas

  use ramses_util
  implicit none

  real(kind=8), parameter :: Msol = 1.9891d+33
  real(kind=8), parameter :: yr = 3.15569d+07
  real(kind=8), parameter :: kparsec = 3.08d+21
  real(kind=8), parameter :: kms = 1.d+05
  real(kind=8), parameter :: pi = 3.14159265d0

  integer :: i,ii,j,ic,choice,np0,n,nb,nstep,nfile,st,rotmat
  integer :: nstar
  integer :: n_frw
  integer,dimension(:),allocatable:: idstar,idj
  real(kind=8) :: aexp,boxlen,scale_l,scale_d,scale_t,scale_v
  real(kind=8) :: scale_cosmict,zz,xc,omega_m,omega_l,omega_k,t,tau
  real(kind=8) :: scale,amass,aexp_ini,time_tot,time_simu,ageu
  real(kind=8),dimension(:,:),allocatable :: xgas,vgas,xstar,vstar,pos,vel
  real(kind=8),dimension(:),allocatable :: zstar,astar,astar2,mstar,mgas,zgas,rhop,pp,tp,dl,mass,met,rho,u,dcell,temp
  real(kind=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  character(len=256) :: arg,repository,fcell,foutgas,fname,fcool,fstar,frot
#ifdef CHEMO
  integer :: nelt
  real(kind=8),dimension(:,:),allocatable :: eltgas,elt
#endif
  integer :: nleaf,lmax,nvarh,iimax
  real(kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax      
  real(kind=8) :: Js,Jxs,Jys,Jzs,mtot,jjx,jjy,jjz,limit,size, &
       ax,bx,cx,ay,by,cy,az,bz,cz,txs,tys,tzs, &
       dx,dy,dd,costh,sinth,cosph,sinph,d1,d2,r0,xc2, &
       theta,phi,theta1,theta2,phi1,phi2,eps
  real(kind=4) :: rji
  integer :: iargc
  character(len=4) :: opt
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

  !! rotmat, output coefficients of the rotation matrix
  !! 0 --> no output 
  !! 1 --> output 
  !! 2 --> output but without rewriting align_gas.dat
  rotmat = 1

  n = iargc()
  if (n < 2) then
     print *, 'usage: align_gas  -nst nstep'
     print *, '                 [-sub subvol_dir] '
     print *, '                 [-out output_dir]'
     print *, '                 [-rji rlim_angmom] '
     print *, '                 [-rot rotmat] '
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
     case ('-nst')
        read(arg,*) nfile
     case ('-out')
        outdir = trim(arg)
     case ('-sub')
        subdir = trim(arg)
     case ('-rji')
        read(arg,*) rji
     case ('-rot')
        read(arg,*) rotmat
     case ('-ver')
        verbose=.true.
     case default
        print '("unknown option ",a4," ignored")', opt
     end select
  end do

  !!======================================================================================

  write(arg,'("/subvol_gas_00",I3,".dat")') nfile
  if (nfile<100) write(arg,'("/subvol_gas_000",I2,".dat")') nfile
  if (nfile<10) write(arg,'("/subvol_gas_0000",I1,".dat")') nfile
  fcell = trim(subdir)//trim(arg)
  print *,'read file: ',trim(fcell)
  open(unit=10,file=fcell,form='unformatted')
  read(10) r0,scale_l,scale_d,scale_t,aexp
  read(10) nleaf
  print *,';leaf =',nleaf
  allocate(xgas(nleaf,3),vgas(nleaf,3),mgas(nleaf),zgas(nleaf),rhop(nleaf),pp(nleaf),dl(nleaf),tp(nleaf))
  read(10) xgas(:nleaf,:3)
  read(10) vgas(:nleaf,:3)
  read(10) mgas(:nleaf)
  read(10) zgas(:nleaf)
  read(10) rhop(:nleaf)
  read(10) pp(:nleaf)
  read(10) dl(:nleaf)
  read(10) tp(:nleaf)
#ifdef CHEMO
  read(10) nelt                                 
  allocate(eltgas(nleaf,nelt))
  read(10) eltgas(:nleaf,:nelt)
#endif
  close(10) 

  if (choice.ne.0) then  

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

  !!======================================================================================
  amass=scale_d*scale_l**3
  amass=amass/Msol
  print '(a8,1pe12.5)','amass =',amass

  size=scale_l/3.08d21
  print *,'size (kpc) =',size

  scale_v=scale_l/scale_t/kms      

  !!======================================================================================
  !! Leo: now align cells

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

  if(rotmat>0)then
     theta1=acos(cz)
     theta2=asin(-cx)
     phi1=acos(by)
     phi2=asin(-ay)
  
     eps=1.e-5
     theta=0.
     if(abs(theta1-theta2)>eps)then
        if(abs(cos(theta2)-cz)<eps)theta=theta2
        if(abs(sin(theta1)+cx)<eps)theta=theta1
     endif
     if(abs(theta1-theta2)<eps)theta=theta1
     if(theta==0.) stop
     
     phi=0.
     if(abs(phi1-phi2)>eps)then
        if(abs(cos(phi2)-by)<eps)phi=phi2 
        if(abs(sin(phi1)+ay)<eps)phi=phi1
     endif
     if(abs(phi1-phi2)<eps)phi=phi1
     if(phi==0.)phi=min(-phi1,pi-phi2)
 
     if(verbose)then 
        print '(1(2x,i5),9(2x,f8.3))',nfile,ax,bx,cx,ay,by,cy,az,bz,cz
        print '(1(2x,i5),6(2x,f8.3))',nfile,theta1,theta2,phi1,phi2
        print '(1(2x,i5),6(2x,f8.3))',nfile,-sin(theta1),cx,cos(theta2),cz
        print *
        print '(1(2x,i5),6(2x,f8.3))',nfile,cos(phi1),by,sin(phi2),-ay
        print '(1(2x,i5),6(2x,f8.3))',nfile,-sin(phi1),ay,cos(phi2),by
        print *
        print '(1(2x,i5),6(2x,f8.3))',nfile,theta,phi,theta*180./pi,phi*180./pi
     endif

     write(arg,'("/rot_mat_00",I3,".dat")') nfile
     if (nfile<100) write(arg,'("/rot_mat_000",I2,".dat")') nfile
     if (nfile<10) write(arg,'("/rot_mat_0000",I1,".dat")') nfile
     frot = trim(subdir)//trim(arg)
     print *,'write file: ',trim(frot)
     open(3, file=frot, status='unknown', form='formatted') 
     write(3,*) 'nfile,ax,bx,cx,ay,by,cy,az,bz,cz,theta(rad),phi(rad),theta(deg),phi(deg)'  
     write(3,'(2x,i5,13(2x,f8.3))') nfile,ax,bx,cx,ay,by,cy,az,bz,cz,theta,phi,theta*180./pi,phi*180./pi
     close(3)
  endif  !if(rotmat)


  if(rotmat==2)then
     deallocate(xgas,vgas,mgas,zgas,rhop,pp,dl,tp)
#ifdef CHEMO
     deallocate(eltgas)
#endif
  endif


  if(rotmat/=2)then
     print *,'ax,bx,cx'
     print *, ax,bx,cx
     print *, ay,by,cy
     print *, az,bz,cz
     
     do j=1,nleaf
        txs=xgas(j,1)
        tys=xgas(j,2)
        tzs=xgas(j,3)
        xgas(j,1)=(ax*txs+bx*tys+cx*tzs)
        xgas(j,2)=(ay*txs+by*tys+cy*tzs)
        xgas(j,3)=(az*txs+bz*tys+cz*tzs)
        
        txs=vgas(j,1)
        tys=vgas(j,2)
        tzs=vgas(j,3)
        vgas(j,1)=(ax*txs+bx*tys+cx*tzs)
        vgas(j,2)=(ay*txs+by*tys+cy*tzs)
        vgas(j,3)=(az*txs+bz*tys+cz*tzs)         
     enddo


     !!======================================================================================
     !! Leo: select a sphere and dump align file

     ! r0 read in subvol header
     allocate(idj(:nleaf))
     ii=0
     idj=0
     do j=1,nleaf
        dd=xgas(j,1)*xgas(j,1)+xgas(j,2)*xgas(j,2)+xgas(j,3)*xgas(j,3)
        if (sqrt(dd)>r0) cycle
        ii=ii+1           
        idj(ii) = j
     enddo
     iimax=ii
     print *,'iimax =',iimax
     
     allocate(pos(:iimax,:3),vel(:iimax,:3),mass(:iimax),met(:iimax),rho(:iimax),u(:iimax),dcell(:iimax),temp(:iimax))
#ifdef CHEMO
     allocate(elt(:iimax,:nelt))
#endif
     do ii=1,iimax
        j = idj(ii)
        pos(ii,1)=xgas(j,1)
        pos(ii,2)=xgas(j,2)
        pos(ii,3)=xgas(j,3)
        vel(ii,1)=vgas(j,1)
        vel(ii,2)=vgas(j,2)
        vel(ii,3)=vgas(j,3)
        mass(ii)=mgas(j)
        met(ii)=zgas(j)
        rho(ii)=rhop(j)
        u(ii)=pp(j)
        dcell(ii)=dl(j)
        temp(ii)=tp(j)
#ifdef CHEMO     
        elt(ii,:nelt)=eltgas(j,:nelt)
#endif
     enddo
     deallocate(xgas,vgas,mgas,zgas,rhop,pp,dl,tp)
     deallocate(idj)
#ifdef CHEMO
     deallocate(eltgas)
#endif
     
     nleaf=iimax
     write(arg,'("/align_gas_00",I3,".dat")') nfile
     if (nfile<100) write(arg,'("/align_gas_000",I2,".dat")') nfile
     if (nfile<10) write(arg,'("/align_gas_0000",I1,".dat")') nfile
     foutgas = trim(outdir)//trim(arg)
     
     open(13, file=foutgas, status='unknown', form='unformatted') 
     rewind 13 
     write(13) r0,scale_l,scale_d,scale_t,aexp
     write(13) nleaf
     write(13) pos(:nleaf,:3)
     write(13) vel(:nleaf,:3)
     write(13) mass(:nleaf)
     write(13) met(:nleaf)
     write(13) rho(:nleaf)
     write(13) u(:nleaf)
     write(13) dcell(:nleaf)
     write(13) temp(:nleaf)
#ifdef CHEMO
     write(13) nelt
     write(13) elt(:nleaf,:nelt)
#endif
     close(13) 
     deallocate(pos,vel,mass,met,rho,u,dcell,temp)
#ifdef CHEMO
     deallocate(elt)
#endif
     print *,'=> dump file: ',trim(foutgas)
     print *
  endif  !if(rotmat/=2)

end program align_gas
