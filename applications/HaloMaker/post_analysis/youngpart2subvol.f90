program main
  !--------------------------------------------------------------------------
  ! from amr2cell.f90:
  ! Ce programme calcule le cube cartesien pour les
  ! variables hydro d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  ! 
  ! and S. Courty's progs
  !
  ! New version by L. Michel-Dansac (feb. 2012)
  ! ce programme lit la position de galaxies/halos et extrait 
  ! les variables hydro d'un sous-volume (cubique) centre sur l'objet
  !
  !--------------------------------------------------------------------------
  implicit none
  integer::ndim,n,i,j,k,twotondim
  integer::ivar,nvar,ncpu,lmax=0
  integer::nx,ny,nz,ilevel,idim,nlevelmax
  integer::ind,ngrida,ngridmax,icpu,ncpu_read
  integer::ihx,ihy
  real(kind=8)::t,aexp,reds,omega_m,omega_l,omega_k,omega_b,boxlen
  integer::imin,imax,jmin,jmax,kmin,kmax
!  integer::nvarh=6
  integer::nvarh
  integer::nx_full,ny_full,nz_full,lmin,nboundary,ngrid_current
  integer::levelmin,nst1,nst2,nb
  integer::ix,iy,iz,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::dx,dx2
  real(KIND=8),dimension(:,:),allocatable::x,xg
  real(KIND=8),dimension(:,:,:),allocatable::var
  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=256)::nomfich,repository,outdir,snapdir,fname,arg,fcool,foutgas
  logical::ok,ok_cell
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  real(KIND=8),dimension(:),allocatable::tp,toto
  real(kind=8)::unit_l,unit_d,unit_t
!#ifdef CHEMO
!  integer:: nelt,ielt
!  real(kind=8),dimension(:,:),allocatable::elt
!#endif

  type level
     integer::ilevel
     integer::ngrid
     !real(KIND=4),dimension(:,:,:),pointer::cube
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
     integer::kmin
     integer::kmax
  end type level
  type(level),dimension(1:100)::grid

  type amr
     real(KIND=8),dimension(:),pointer::dcell
     real(KIND=8),dimension(:,:),pointer::hydro,pos
  end type amr
  type(amr) :: cell

  type cooling_table
     real(kind=8),dimension(:)    ,pointer::nH
     real(kind=8),dimension(:)    ,pointer::T2
     real(kind=8),dimension(:)    ,pointer::T2eq
     real(kind=8),dimension(:,:)  ,pointer::metal
     real(kind=8),dimension(:,:)  ,pointer::cool
     real(kind=8),dimension(:,:)  ,pointer::heat
     real(kind=8),dimension(:,:)  ,pointer::metal_prime
     real(kind=8),dimension(:,:)  ,pointer::cool_prime
     real(kind=8),dimension(:,:)  ,pointer::heat_prime
     real(kind=8),dimension(:,:)  ,pointer::mu
     real(kind=8),dimension(:,:,:),pointer::spec
  end type cooling_table
  type(cooling_table)::cooling

  real(kind=8)::wxa,wxb,wya,wyb,muval,xx,pre,rho,dyy,yy,dxx
  integer :: nleaf,nleafcpu,ic,nfile,nstep,n1,n2,ii
  real(kind=8)::size,H0,xlim,xlimb,rvir,mvir
  real(kind=8),dimension(3)::cm
  integer, dimension(:),allocatable :: tabnfile
  real(kind=8) :: x0,y0,z0
  real(kind=8), dimension(:,:),allocatable :: poscm

  real(kind=8) :: r0,scale_l,scale_d,scale_t,scale_v,ageu,time_simu,dzz,dr_cell
  integer :: istar,np,ist,nstar
  integer, dimension(:),allocatable :: itmp
  real(kind=8),dimension(:,:),allocatable :: xstar,vstar,xgstar
  real(kind=8),dimension(:),allocatable   :: mgstar,zgstar,astar,astar2
  character(len=256) :: fstar,subdir


  !! set parameters by default
  outdir = "."
  snapdir = "."
  nst1 = 1
  nst2 = 0
  xlim = 200.  !! in kpc
  subdir = "."

  call read_params

  !-------------------------------------------
  ! read new center of mass of halos/galaxies
  !-------------------------------------------
  open(13,file=fname,form='formatted')
  read(13,*) nb
  print *,'nb ',nb
  allocate(tabnfile(nb),poscm(nb,3))
  do n=1,nb
     read(13,*) j,nfile,nstep,aexp,reds,ic,x0,y0,z0,rvir,mvir
     print '(3(2x,i5),2x,f7.4,2x,f6.3,2x,i12,4(2x,f9.6),2x,e11.4)',j,nfile,nstep,aexp,reds,ic,x0,y0,z0,rvir,mvir
     tabnfile(n)=nfile
     poscm(n,1)=x0
     poscm(n,2)=y0
     poscm(n,3)=z0
  enddo
  close(13)
  if (nst2 == 0) nst2=nb
  print *

  !---------------------
  ! Loop over snapshots
  !---------------------
  do n=nst1,nst2
     nfile=tabnfile(n)
     cm(1)=poscm(n,1)
     cm(2)=poscm(n,2)
     cm(3)=poscm(n,3)         
     print *,'n ',n,nfile
     print *,'x0,y0,z0 ',cm(1),cm(2),cm(3)

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
  allocate(xstar(nstar,3),vstar(nstar,3),astar2(nstar))
  read(13) xstar(:nstar,:3)
  read(13) vstar(:nstar,:3)
  read(13) astar2(:nstar)
  read(13) astar2(:nstar)
  read(13) astar2(:nstar)
  close(13)

print *,ageu,time_simu,maxval(astar2(:nstar))

  xstar(:nstar,1)=xstar(:nstar,1)+0.5
  xstar(:nstar,2)=xstar(:nstar,2)+0.5
  xstar(:nstar,3)=xstar(:nstar,3)+0.5
  do j=1,3
     dx=cm(j)-0.5
     xstar(:nstar,j)=xstar(:nstar,j)+dx 
     where (xstar(:nstar,j) < 0.) xstar(:nstar,j)=xstar(:nstar,j)+1.
     where (xstar(:nstar,j) > 1.) xstar(:nstar,j)=xstar(:nstar,j)-1.
  enddo

  ii=0
  allocate(itmp(nstar))
  itmp=0
  do j=1,nstar
!???
     if(astar2(j)<(maxval(astar2(:nstar))-5.e8/1.e9))cycle
     ii=ii+1
     itmp(ii)=j
  enddo
  np=ii
  print *,'nb of young stars ',np
  allocate(xgstar(np,3),mgstar(np),zgstar(np))
  !!======================================================================================

     write(nchar,'("00",I3)') nfile
     if (nfile<100) write(nchar,'("000",I2)') nfile
     if (nfile<10) write(nchar,'("0000",I1)') nfile
     arg='/output_'//trim(nchar)
     repository = trim(snapdir)//trim(arg)
     print *,'output directory = ',trim(repository)

     !-----------------------------------------------
     ! Lecture du fichier hydro au format RAMSES
     !-----------------------------------------------
     !ipos=INDEX(repository,'output_')
     !nchar=repository(ipos+7:ipos+13)
     nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
     inquire(file=nomfich, exist=ok) ! verify input file 
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found.'
        stop
     endif
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
     inquire(file=nomfich, exist=ok) ! verify input file 
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found.'
        stop
     endif

     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
     open(unit=10,file=nomfich,status='old',form='unformatted')
     read(10)ncpu
     read(10)ndim
     read(10)nx,ny,nz
     read(10)nlevelmax
     read(10)ngridmax
     read(10)nboundary
     read(10)ngrid_current
     read(10)boxlen
     close(10)
     twotondim=2**ndim
     xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

     allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
     allocate(ngridlevel(1:ncpu,1:nlevelmax))
     if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

     if(ndim==2)then
        write(*,*)'Output file contains 2D data'
        write(*,*)'Aborting'
        stop
     endif

     nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
     inquire(file=nomfich, exist=ok) ! verify input file
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found.'
        stop
     endif
     open(unit=10,file=nomfich,form='formatted',status='old')
     read(10,*)
     read(10,*)
     read(10,'("levelmin    =",I11)')levelmin
     read(10,*)
     read(10,*)
     read(10,*)
     read(10,*)

     read(10,*)
     read(10,'("time        =",E23.15)')t
     read(10,'("aexp        =",E23.15)')aexp
     read(10,'("H0          =",E23.15)')H0
     read(10,'("omega_m     =",E23.15)')omega_m
     read(10,'("omega_l     =",E23.15)')omega_l
     read(10,'("omega_k     =",E23.15)')omega_k
     read(10,'("omega_b     =",E23.15)')omega_b
     read(10,'("unit_l      =",E23.15)')unit_l
     read(10,'("unit_d      =",E23.15)')unit_d
     read(10,'("unit_t      =",E23.15)')unit_t
     read(10,*)

     read(10,'("ordering type=",A80)'),ordering
     write(*,'(" ordering type=",A20)'),TRIM(ordering)
     read(10,*)
     allocate(cpu_list(1:ncpu))
     if(TRIM(ordering).eq.'hilbert')then
        allocate(bound_key(0:ncpu))
        allocate(cpu_read(1:ncpu))
        cpu_read=.false.
        do impi=1,ncpu
           read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
        end do
     endif
     close(10)

     print *,'cosmology ',omega_m,omega_l,omega_k,omega_b,H0
     print *,unit_l,unit_d,unit_t
     size=unit_l/3.08d21
     print *,'size (kpc) ',size


     !--------------------------------------------
     ! define position and size of the sub-volume
     !--------------------------------------------
     size=dble(unit_l)/3.08d21
     !!print *,'size (kpc) ',size
     print *,'sub-volume position and size'
     print *,'cm =',cm(1),cm(2),cm(3)
     xlimb = xlim/size
     print '(a7,2x,e11.4,2x,f6.2)','xlim =',xlimb,xlim
     xmin = cm(1)-xlimb
     xmax = cm(1)+xlimb
     ymin = cm(2)-xlimb
     ymax = cm(2)+xlimb
     zmin = cm(3)-xlimb
     zmax = cm(3)+xlimb
     print *,xmin,xmax
     print *,ymin,ymax
     print *,zmin,zmax

     if(lmax==0)then
        lmax=nlevelmax
     endif
     write(*,*)'time=',t
     write(*,*)'Working resolution =',2**lmax

     if(TRIM(ordering).eq.'hilbert')then

        dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
        do ilevel=1,lmax
           dx=0.5d0**ilevel
           if(dx.lt.dmax)exit
        end do
        lmin=ilevel
        bit_length=lmin-1
        maxdom=2**bit_length
        imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
        if(bit_length>0)then
           imin=int(xmin*dble(maxdom))
           imax=imin+1
           jmin=int(ymin*dble(maxdom))
           jmax=jmin+1
           kmin=int(zmin*dble(maxdom))
           kmax=kmin+1
        endif

        dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
        ndom=1
        if(bit_length>0)ndom=8
        idom(1)=imin; idom(2)=imax
        idom(3)=imin; idom(4)=imax
        idom(5)=imin; idom(6)=imax
        idom(7)=imin; idom(8)=imax
        jdom(1)=jmin; jdom(2)=jmin
        jdom(3)=jmax; jdom(4)=jmax
        jdom(5)=jmin; jdom(6)=jmin
        jdom(7)=jmax; jdom(8)=jmax
        kdom(1)=kmin; kdom(2)=kmin
        kdom(3)=kmin; kdom(4)=kmin
        kdom(5)=kmax; kdom(6)=kmax
        kdom(7)=kmax; kdom(8)=kmax

        do i=1,ndom
           if(bit_length>0)then
              call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           else
              order_min=0.0d0
           endif
           bounding_min(i)=(order_min)*dkey
           bounding_max(i)=(order_min+1.0D0)*dkey
        end do

        cpu_min=0; cpu_max=0
        do impi=1,ncpu
           do i=1,ndom
              if (   bound_key(impi-1).le.bounding_min(i).and.&
                   & bound_key(impi  ).gt.bounding_min(i))then
                 cpu_min(i)=impi
              endif
              if (   bound_key(impi-1).lt.bounding_max(i).and.&
                   & bound_key(impi  ).ge.bounding_max(i))then
                 cpu_max(i)=impi
              endif
           end do
        end do

        ncpu_read=0
        do i=1,ndom
           do j=cpu_min(i),cpu_max(i)
              if(.not. cpu_read(j))then
                 ncpu_read=ncpu_read+1
                 cpu_list(ncpu_read)=j
                 cpu_read(j)=.true.
              endif
           enddo
        enddo
     else
        ncpu_read=ncpu
        do j=1,ncpu
           cpu_list(j)=j
        end do
     end  if

     !-----------------------------
     ! Compute hierarchy
     !-----------------------------
     do ilevel=1,lmax
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel
        imin=int(xmin*dble(nx_full))+1
        imax=int(xmax*dble(nx_full))+1
        jmin=int(ymin*dble(ny_full))+1
        jmax=int(ymax*dble(ny_full))+1
        kmin=int(zmin*dble(nz_full))+1
        kmax=int(zmax*dble(nz_full))+1
        grid(ilevel)%imin=imin
        grid(ilevel)%imax=imax
        grid(ilevel)%jmin=jmin
        grid(ilevel)%jmax=jmax
        grid(ilevel)%kmin=kmin
        grid(ilevel)%kmax=kmax
     end do

     print *,'levelmin ',levelmin
     print *,'levelmax ',lmax
     print *,'ncpu_read =',ncpu_read,ncpu

     ! Loop over processor files
     do k=1,ncpu_read
        icpu=cpu_list(k)
        call title(icpu,ncharcpu)

        ! Open AMR file and skip header
        nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
        open(unit=10,file=nomfich,status='old',form='unformatted')
        !write(*,*)'Processing file '//TRIM(nomfich)
        do i=1,21
           read(10)
        end do
        ! Read grid numbers
        read(10)ngridlevel
        ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
        read(10)
        if(nboundary>0)then
           do i=1,2
              read(10)
           end do
           read(10)ngridbound
           ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
        endif
        read(10)
        ! ROM: comment the single follwing line for old stuff
        read(10)
        if(TRIM(ordering).eq.'bisection')then
           do i=1,5
              read(10)
           end do
        else
           read(10)
        endif
        read(10)
        read(10)
        read(10)

        ! Open HYDRO file and skip header
        nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
        open(unit=11,file=nomfich,status='old',form='unformatted')
        read(11)
        read(11)nvarh
        read(11)
        read(11)
        read(11)
        read(11)

        ! Loop over levels
        do ilevel=1,lmax

           ! Geometry
           dx=0.5**ilevel
           dx2=0.5*dx
!           nx_loc=(icoarse_max-icoarse_min+1)
!           scale=boxlen/dble(nx_loc)
!           dx_loc=dx*scale
           nx_full=2**ilevel
           ny_full=2**ilevel
           nz_full=2**ilevel
           do ind=1,twotondim
              iz=(ind-1)/4
              iy=(ind-1-4*iz)/2
              ix=(ind-1-2*iy-4*iz)
              xc(ind,1)=(dble(ix)-0.5D0)*dx
              xc(ind,2)=(dble(iy)-0.5D0)*dx
              xc(ind,3)=(dble(iz)-0.5D0)*dx
           end do

           ! Allocate work arrays
           ngrida=ngridfile(icpu,ilevel)
           grid(ilevel)%ngrid=ngrida
           if(ngrida>0)then
              allocate(xg(1:ngrida,1:ndim))
              allocate(son(1:ngrida,1:twotondim))
              allocate(var(1:ngrida,1:twotondim,1:nvarh))
              allocate(x  (1:ngrida,1:ndim))
              !allocate(rho(1:ngrida))
              allocate(ref(1:ngrida))
           endif

           ! Loop over domains
           do j=1,nboundary+ncpu

              ! Read AMR data
              if(ngridfile(j,ilevel)>0)then
                 read(10) ! Skip grid index
                 read(10) ! Skip next index
                 read(10) ! Skip prev index
                 ! Read grid center
                 do idim=1,ndim
                    if(j.eq.icpu)then
                       read(10)xg(:,idim)
                    else
                       read(10)
                    endif
                 end do
                 read(10) ! Skip father index
                 do ind=1,2*ndim
                    read(10) ! Skip nbor index
                 end do
                 ! Read son index
                 do ind=1,twotondim
                    if(j.eq.icpu)then
                       read(10)son(:,ind)
                    else
                       read(10)
                    end if
                 end do
                 ! Skip cpu map
                 do ind=1,twotondim
                    read(10)
                 end do
                 ! Skip refinement map
                 do ind=1,twotondim
                    read(10)
                 end do
              endif

              ! Read HYDRO data
              read(11)
              read(11)
              if(ngridfile(j,ilevel)>0)then
                 ! Read hydro variables
                 do ind=1,twotondim
                    do ivar=1,nvarh
                       if(j.eq.icpu)then
                          read(11)var(:,ind,ivar)
                       else
                          read(11)
                       end if
                    end do
                 end do
              end if
           end do

           ! Compute map
           if(ngrida>0)then

              ! Loop over cells
              do ind=1,twotondim

                 ! Compute cell center
                 do i=1,ngrida
                    x(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                    x(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                    x(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
                 end do
                 ! Check if cell is refined
                 do i=1,ngrida
                    ref(i)=son(i,ind)>0.and.ilevel<lmax
                 end do
                 ! Store data cube
                 do i=1,ngrida
                    ok_cell= .not.ref(i).and. &
                         & (x(i,1)+dx2)>=xmin.and.&
                         & (x(i,2)+dx2)>=ymin.and.&
                         & (x(i,3)+dx2)>=zmin.and.&
                         & (x(i,1)-dx2)<=xmax.and.&
                         & (x(i,2)-dx2)<=ymax.and.&
                         & (x(i,3)-dx2)<=zmax
                    if(ok_cell)then

                       do ist=1,np
                          istar=itmp(ist)
                          dxx=x(i,1)-xstar(istar,1)
                          dyy=x(i,2)-xstar(istar,2)
                          dzz=x(i,3)-xstar(istar,3)
                          dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))

print *,ist,astar2(istar),xstar(istar,1),xstar(istar,2),xstar(istar,3)

                          if(dr_cell.le.dx/2.0)then
                             do idim=1,ndim
                                xgstar(istar,idim)=x(i,idim)
                             enddo
                             mgstar(istar)=var(i,ind,1)*dx*dx*dx
                             zgstar(istar)=var(i,ind,6)

print *,ist,xgstar(istar,1),xgstar(istar,2),xgstar(istar,3)
print *,ist,mgstar(istar),zgstar(istar)

                          endif
                       enddo

                    end if
                 end do

              end do
              ! End loop over cell

              deallocate(xg,son,var,ref,x)
           endif

        end do
        ! End loop over levels

        close(10)
        close(11)

     end do
     ! End loop over cpus

     !--------------------------------------------
     ! recenter cells
     !--------------------------------------------
     do j=1,3
        dx=cm(j)-0.5
        xgstar(:np,j)=xgstar(:np,j)-dx
        where (xgstar(:np,j) < 0.) xgstar(:np,j)=xgstar(:np,j)+1.
        where (xgstar(:np,j) > 1.) xgstar(:np,j)=xgstar(:np,j)-1.
     enddo
     do j=1,3
        xgstar(:np,j)=xgstar(:np,j)-0.5
     enddo

     !--------------------------------------------
     ! write file
     !--------------------------------------------

     arg='/subvol_gasstar_'//trim(nchar)//'.dat'
     foutgas = trim(outdir)//trim(arg)

     open(13, file=foutgas, status='unknown', form='unformatted') 
     write(13) r0,scale_l,scale_d,scale_t,aexp,ageu,time_simu
     write(13) np
     write(13) (itmp(i),i=1,np)
     write(13) ((xgstar(i,j),i=1,np),j=1,3)   !position
     write(13) (mgstar(i),i=1,np)             ! mass
     write(13) (zgstar(i),i=1,np)             ! metallicity
     close(13) 

     print *,'=> dump file: ',trim(foutgas)
     print *

     deallocate(cpu_list)
     deallocate(bound_key)
     deallocate(cpu_read)
     deallocate(ngridfile)
     deallocate(ngridlevel)
     if(nboundary>0)deallocate(ngridbound)

     deallocate(xstar,vstar,astar2)
     deallocate(itmp,xgstar,mgstar,zgstar)
  enddo
  ! end loop over snapshots

  deallocate(tabnfile,poscm)


!=======================================================================
!=======================================================================
contains
  
  subroutine read_params
    
    implicit none
    
    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok
    
    n = iargc()

    if (n < 2) then
       print *, 'usage: amr2subvol  -inp input_file'
       print *, '                  [-ns1 nstep_first]'
       print *, '                  [-ns2 nstep_last]'
       print *, '                  [-sna snapshot_dir]'
       print *, '                  [-out output_dir]'
       print *, '                  [-xli xlim] '
       print *, '                  [-sub subvol_dir] '
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
       case ('-inp')
          fname = trim(arg)
       case ('-out')
          outdir = trim(arg)
       case ('-sna')
          snapdir = trim(arg)
       case ('-ns1')
          read(arg,*) nst1
       case ('-ns2')
          read(arg,*) nst2
       case ('-xli')
          read(arg,*) xlim
       case ('-sub')
          subdir = trim(arg)
       case default
          print '("unknown option ",a4," ignored")', opt
       end select
    enddo

    return
    
  end subroutine read_params
  
end program main

!=======================================================================
!=======================================================================
subroutine title(n,nchar)

  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title

!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d
!================================================================
!================================================================
!================================================================
!================================================================
