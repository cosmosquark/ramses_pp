program ccmassgal
  !----------------------------------------------------------------------------------
  ! from part2map.f90:
  ! Ce programme calcule la carte de densite surfacique projetee
  ! des particules de matiere noire d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  !
  ! and cdm_loop.dat (S. Courty)
  ! 
  ! Version par L. Michel-Dansac (feb. 2012)
  ! ce programme recalcule le centre de masse sur les etoiles 
  ! de la galaxie au centre du halo
  ! Usage: ccmassgal -inp progenitor.dat [ -ns1 XX -ns2 XX -sna snapdir -out outdir -rli 5.]
  ! --> galaxies_ccmass.dat
  ! ns1 and ns2 are relative to progenitor.dat, not to snapshot number...
  !-----------------------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,n,i,j,k,icpu,ipos,nstar
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
  integer::idim,ncpu_read
  real(KIND=8)::mtot
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin
  real(KIND=8)::dx,dy,deltax,boxlen
  real(KIND=8)::aexp,t,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d
  real(KIND=8),dimension(:,:),allocatable::x,v
  real(KIND=8),dimension(:)  ,allocatable::m,age,z
  integer,dimension(:)  ,allocatable::id
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering,format_grille
  character(LEN=256)::nomfich,repository,outfich,arg,fname,outfile,snapdir,outdir
  logical::ok,ok_part,periodic=.false.,star=.false.,ok_star
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  integer::npst,npdm,npstcpu,npdmcpu,ii,n1,n2,nfile,nb,nstep,ic
  real(KIND=8)::size,rlim,dd,dx0,dy0,dz0
  real(KIND=8),dimension(:,:),allocatable::xstar,vstar
  real(KIND=8),dimension(:)  ,allocatable::mstar,astar,zstar,astar2
  integer,dimension(:)  ,allocatable::idstar
  real(kind=8),dimension(3)::cm,cm2
  real(kind=8), dimension(:,:),allocatable :: tdiff
  integer, dimension(:),allocatable :: tabnfile
  real(kind=8) :: x0,y0,z0,xc,xc2
  real(kind=8), dimension(:,:),allocatable :: poscm,poscm2
  real(kind=8) :: rvir,mvir,xlim

  !! set parameters by default
  rlim = 3.   !! in kpc
  outdir = "."
  snapdir = "."
  n1 = 1
  n2 = 0
  xlim = 30.  !! in kpc

  call read_params

  !-------------------------------
  ! read center of mass of halos
  !-------------------------------

  open(13,file=fname,form='formatted')
  read(13,*) nb
  print *,'nb ',nb
  allocate(tabnfile(nb),poscm(nb,3))
  do n=1,nb
     read(13,*) j,nfile,nstep,xc,xc2,ic,x0,y0,z0,rvir,mvir
     print '(3(2x,i5),2x,f7.4,2x,f6.3,2x,i12,4(2x,f9.6),2x,e11.4)',j,nfile,nstep,xc,xc2,ic,x0,y0,z0,rvir,mvir
     tabnfile(n)=nfile
     poscm(n,1)=x0
     poscm(n,2)=y0
     poscm(n,3)=z0
  enddo
  close(13)
  allocate(poscm2(nb,3),tdiff(nb,3))
  if (n2 == 0) n2=nb

  !------------------------------------------------------
  ! Loop over snapshots
  !------------------------------------------------------
  do n=n1,n2                     
     nfile=tabnfile(n)
     cm(1)=poscm(n,1)
     cm(2)=poscm(n,2)
     cm(3)=poscm(n,3)         
     print *,'n ',n,nfile
     print *,'x0,y0,z0 ',cm(1),cm(2),cm(3)


     write(nchar,'("00",I3)') nfile
     if (nfile<100) write(nchar,'("000",I2)') nfile
     if (nfile<10) write(nchar,'("0000",I1)') nfile
     arg='/output_'//trim(nchar)
     repository = trim(snapdir)//trim(arg)
     print *,'output directory = ',trim(repository)

     !-----------------------------------------------
     ! Lecture du fichier particules au format RAMSES
     !-----------------------------------------------
     !ipos=INDEX(repository,'output_')
     !nchar=repository(ipos+7:ipos+13)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
     inquire(file=nomfich, exist=ok) ! verify input file 
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found.'
        stop
     endif

     nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
     inquire(file=nomfich, exist=ok) ! verify input file 
     if ( .not. ok ) then
        print *,TRIM(nomfich)//' not found.'
        stop
     endif
     open(unit=10,file=nomfich,form='formatted',status='old')
     read(10,'("ncpu        =",I11)')ncpu
     read(10,'("ndim        =",I11)')ndim
     read(10,'("levelmin    =",I11)')levelmin
     read(10,'("levelmax    =",I11)')levelmax
     read(10,*)
     read(10,*)
     read(10,*)

     read(10,'("boxlen      =",E23.15)')boxlen
     read(10,'("time        =",E23.15)')t
     read(10,'("aexp        =",E23.15)')aexp
     read(10,'("H0          =",E23.15)')h0
     read(10,'("omega_m     =",E23.15)')omega_m
     read(10,'("omega_l     =",E23.15)')omega_l
     read(10,'("omega_k     =",E23.15)')omega_k
     read(10,'("omega_b     =",E23.15)')omega_b
     read(10,'("unit_l      =",E23.15)')unit_l
     read(10,'("unit_d      =",E23.15)')unit_d
     read(10,'("unit_t      =",E23.15)')unit_t
     read(10,*)

     read(10,'("ordering type=",A80)') ordering
     write(*,'(" ordering type = ",A20)') TRIM(ordering)
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

     !--------------------------------------------
     ! define position and size of the sub-volume
     !--------------------------------------------
     size=dble(unit_l)/3.08d21
     !!print *,'size (kpc) ',size
     print *,'sub-volume position and size'
     print *,'cm =',cm(1),cm(2),cm(3)
     dx = xlim/size
     print *,'dx =',dx,dx*size
     xmin = cm(1)-dx
     xmax = cm(1)+dx
     ymin = cm(2)-dx
     ymax = cm(2)+dx
     zmin = cm(3)-dx
     zmax = cm(3)+dx
     print *,xmin,xmax
     print *,ymin,ymax
     print *,zmin,zmax


     if(TRIM(ordering).eq.'hilbert')then

        dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
        do ilevel=1,levelmax
           deltax=0.5d0**ilevel
           if(deltax.lt.dmax)exit
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

        dkey=(dble(2**(levelmax+1)/dble(maxdom)))**ndim
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

     print *,'ncpu_read =',ncpu_read,ncpu
     print *,'boxlen =',boxlen

!!$     npart=0
!!$     do k=1,ncpu_read
!!$        icpu=cpu_list(k)
!!$        call title(icpu,ncharcpu)
!!$        nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
!!$        open(unit=1,file=nomfich,status='old',form='unformatted')
!!$        read(1)ncpu2
!!$        read(1)ndim2
!!$        read(1)npart2
!!$        read(1)
!!$        read(1)nstar
!!$        close(1)
!!$        npart=npart+npart2
!!$     end do
!!$     write(*,*)'Found ',npart,' particles.'
!!$     if(nstar>0)then
!!$        if(star)then
!!$           write(*,*)'Keeping star particles.'
!!$        else
!!$           write(*,*)'Discard star particles.'
!!$        endif
!!$     endif

     ! find nstar inside sub-volume
     npdm=0
     npst=0
     npart=0
     do k=1,ncpu_read
        icpu=cpu_list(k)
        call title(icpu,ncharcpu)
        nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
        open(unit=1,file=nomfich,status='old',form='unformatted')
        !     write(*,*)'Processing file '//TRIM(nomfich)
        read(1)ncpu2
        read(1)ndim2
        read(1)npart2
        read(1)
        read(1)nstar
        npart=npart+npart2
        !! Leo: nstar is the total number of star whereas npart2 is the number of part (ndm_cpu+nstar_cpu) in this cpu file...
        read(1)
        read(1)
        read(1)
        allocate(m(1:npart2))
        if(nstar>0)then
           allocate(age(1:npart2))
           allocate(id(1:npart2))
        endif
        allocate(x(1:npart2,1:ndim2))
        ! Read position
        do i=1,ndim
           read(1)m
           x(1:npart2,i)=m/boxlen
        end do
        ! Skip velocity
        do i=1,ndim
           read(1)m
        end do
        ! Read mass
        read(1)m
        if(nstar>0)then
           read(1)id
           read(1) ! Skip level
           read(1)age
        endif
        close(1)

        !! count particles inside sub-volume
        npdmcpu=0
        npstcpu=0
        do i=1,npart2
           ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
                &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
                &   x(i,3)>=zmin.and.x(i,3)<=zmax)
           if(nstar>0)then
              ! Leo:
              ! Ramses_v3.07 public version: id=0 are GMCs locked mass -> do not select id=0
              !ok_star=ok_part.and.(age(i).ne.0.0d0).and.(id(i).gt.0)
              ! Ramses v3.07 new version: very young stars do have id<0 -> select id<0
              ok_star=ok_part.and.(age(i).ne.0.0d0)           
              ok_part=ok_part.and.(age(i)==0.0d0)
           endif
           if(ok_part)then
              npdm=npdm+1
              npdmcpu=npdmcpu+1
           endif
           if(ok_star)then
              npst=npst+1
              npstcpu=npstcpu+1
           endif
        enddo
        !!print *,npdmcpu,npstcpu
        deallocate(x,m)
        if(nstar>0)deallocate(age,id)
     end do
     write(*,*)'Found ',npart,' particles.'
     print *,'Ndm_subvol =',npdm
     print *,'Nst_subvol =',npst

     !! read position & mass for stars inside sub-volume
     !!allocate(xstar(npst,3),vstar(npst,3),mstar(npst),zstar(npst),astar(npst),astar2(npst),idstar(npst))
     allocate(xstar(npst,3),mstar(npst))
     npdmcpu=0
     npstcpu=0
     do k=1,ncpu_read
        icpu=cpu_list(k)
        call title(icpu,ncharcpu)
        nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
        open(unit=1,file=nomfich,status='old',form='unformatted')
        !     write(*,*)'Processing file '//TRIM(nomfich)
        read(1)ncpu2
        read(1)ndim2
        read(1)npart2
        read(1)
        read(1)nstar
        !! Leo: nstar is the total number of star whereas npart2 is the number of part (dm_cpu+nstar_cpu) in this cpu file...
        read(1)
        read(1)
        read(1)
        allocate(m(1:npart2))
        allocate(id(1:npart2))
        if(nstar>0)then
           allocate(age(1:npart2))
        endif
        allocate(x(1:npart2,1:ndim2))
        ! Read position
        do i=1,ndim
           read(1)m
           x(1:npart2,i)=m/boxlen
        end do
        ! Skip velocity
        do i=1,ndim
           read(1)m
        end do
        ! Read mass
        read(1)m
        ! Read identity
        read(1)id
        ! Skip level
        read(1) ! skip level
        if(nstar>0)then
           ! Read birth epoch 
           read(1)age
        endif
        close(1)
        do i=1,npart2
           ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
                &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
                &   x(i,3)>=zmin.and.x(i,3)<=zmax)
           if(nstar>0)then
              ! Leo:
              ! Ramses_v3.07 public version: id=0 are GMCs locked mass -> do not select id=0
              !ok_star=ok_part.and.(age(i).ne.0.0d0).and.(id(i).gt.0)
              ! Ramses v3.07 new version: very young stars do have id<0 -> select id<0
              ok_star=ok_part.and.(age(i).ne.0.0d0)           
              ok_part=ok_part.and.(age(i)==0.0d0)
           endif
           if(ok_part)then
              npdmcpu=npdmcpu+1
           endif
           if(ok_star)then
              npstcpu=npstcpu+1
              do j=1,ndim2
                 xstar(npstcpu,j)=x(i,j)*boxlen
              enddo
              mstar(npstcpu)=m(i)
           endif
        enddo
        !!print *,npdmcpu,npstcpu
        deallocate(x,m,id)
        if(nstar>0)deallocate(age)
     end do
     print *,npdmcpu,npstcpu

     !---------------------------------
     ! check CoM of galaxy using stars
     !---------------------------------

     !sphere radius, kpc
     print *,'rlim =',rlim,rlim/size

     !New Center of Mass
     cm2=0.d0
     ii=0
     mtot=0.d0

     !  print *,minval(abs(xstar(:,1)-cm(1))),maxval(abs(xstar(:,1)-cm(1)))
     !  print *,minval(abs(xstar(:,2)-cm(2))),maxval(abs(xstar(:,2)-cm(3)))
     !  print *,minval(abs(xstar(:,3)-cm(3))),maxval(abs(xstar(:,2)-cm(3)))
     !  print *,minval(xstar(:,2)),maxval(xstar(:,2)),cm(2)
     !  print *,minval(xstar(:,3)),maxval(xstar(:,3)),cm(3)

     do j=1,npst
        !!if (idstar(j)==0) cycle
        dx0=cm(1)-xstar(j,1)
        dy0=cm(2)-xstar(j,2)
        dz0=cm(3)-xstar(j,3)
        dd=sqrt(dx0*dx0+dy0*dy0+dz0*dz0)

        if (dd>(rlim/size)) cycle

        ii=ii+1
        mtot=mtot+mstar(j)
        do i=1,3 
           cm2(i)=cm2(i)+xstar(j,i)*mstar(j)
        enddo
     enddo
     print *,'ii ',ii,mtot
     if (ii == 0) then
        print *,'no stars in this galaxy...'
        stop
     endif
     cm2(:3)=cm2(:3)/mtot
     print *,'cm2 ',cm2(1),cm2(2),cm2(3)
     print *,'diff ',abs(cm(1)-cm2(1))*size,abs(cm(2)-cm2(2))*size,abs(cm(3)-cm2(3))*size
     print *

     tdiff(n,1)=abs(cm(1)-cm2(1))*size
     tdiff(n,2)=abs(cm(2)-cm2(2))*size
     tdiff(n,3)=abs(cm(3)-cm2(3))*size

     poscm2(n,1)=cm2(1)
     poscm2(n,2)=cm2(2)
     poscm2(n,3)=cm2(3)

     deallocate(xstar,mstar)
     deallocate(cpu_list)
     deallocate(bound_key)
     deallocate(cpu_read)


  enddo

  do n=n1,n2
     print '(i12,3(2x,f11.8))',n,poscm(n,1),poscm(n,2),poscm(n,3)
     print '(12x,3(2x,f11.8))',poscm2(n,1),poscm2(n,2),poscm2(n,3)
     print '(12x,3(2x,e11.4))',tdiff(n,1),tdiff(n,2),tdiff(n,3)
  enddo

  open(13,file=fname,form='formatted')
  read(13,*) nb

  outfile=trim(outdir)//'/galaxies_ccmass.dat'
  open(10,file=outfile,form='formatted')
  write(10,*) n2
      
  do n=1,n2
     read(13,*) j,nfile,nstep,xc,xc2,ic,x0,y0,z0,rvir,mvir
     if (n>=n1) then 
        write(10,'(3(2x,i5),2x,f7.4,2x,f6.3,2x,i12,4(2x,f10.7),2x,e11.4)') j,nfile,nstep,xc,xc2,ic, &
             poscm2(n,1),poscm2(n,2),poscm2(n,3),rvir,mvir
     else 
        write(10,'(3(2x,i5),2x,f7.4,2x,f6.3,2x,i12,4(2x,f10.7),2x,e11.4)') j,nfile,nstep,xc,xc2,ic,&
             x0,y0,z0,rvir,mvir
     endif
  enddo
  close(10)
  close(13)
      
  deallocate(tabnfile,poscm,poscm2)


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
       print *, 'usage: ccmassgal  -inp input_file'
       print *, '                [-ns1 nstep_first]'
       print *, '                [-ns2 nstep_last]'
       print *, '                [-sna snapshot_dir]'
       print *, '                [-out output_dir]'
       print *, '                [-rli rlim_com] '
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
          read(arg,*) n1
       case ('-ns2')
          read(arg,*) n2
       case ('-rli')
          read(arg,*) rlim
       case default
          print '("unknown option ",a4," ignored")', opt
       end select
    enddo

    return

  end subroutine read_params

end program ccmassgal

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
subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min,age_tot
  real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

!  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
!     write(*,*)'Error: non-physical cosmological constants'
!     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
!     write(*,*)'The sum must be equal to 1.0, but '
!     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
!     stop
!  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0
  
  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau
     
     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
  end do

  age_tot=-t
  write(*,666)-t
  666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable
  
  axp_t = 1.d0
  t = 0.d0
  axp_tau = 1.d0
  tau = 0.d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if

  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman
!=======================================================================
!=======================================================================
function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau
!=======================================================================
!=======================================================================
function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt
!=======================================================================
!=======================================================================
!=======================================================================
