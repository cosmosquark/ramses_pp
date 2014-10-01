!--------------------------------------------------------------------------------------------------------------
! programme treebrick2ascii.f90
! Purpose: read treebrick file from HaloMaker and write two ascii files
! ~/Ramses_v3.01/Stat/read_treebrick -inp output_00030 halomaker_00030.nod halomaker_00030.nei tree_bricks030 
! from S. Courty original version
! Version by L. Michel-Dansac (feb.2012)
!--------------------------------------------------------------------------------------------------------------
module tree_defs

  type vector
     real(kind=4)            :: x,y,z
  end type vector

  type father
     real(kind=4),pointer    :: mass_fathers(:)  ! percentage of father mass that goes into son halo
     integer(kind=4),pointer :: list_fathers(:)
!     integer(kind=4) :: list_fathers(:)
     integer(kind=4)         :: nb_fathers
  end type father        
                             
  type shape                 
     real(kind=4)            :: a,b,c
  end type shape             
                             
  type baryon                
     real(kind=4)            :: rvir,mvir,tvir,cvel
  end type baryon

  type son
     integer(kind=4),pointer :: list_sons(:)
     integer(kind=4)         :: nb_sons
  end type son               
                             
  type hprofile              
     real(kind=4)            :: rho_0,r_c
  end type hprofile          
                             
  type halo                  
     type (baryon)           :: datas
     type (father)           :: my_fathers
     type (son)              :: my_sons
     type (shape)            :: sh
     type (vector)           :: p
     type (vector)           :: v
     type (vector)           :: L
     type (hprofile)         :: halo_profile 
     integer(kind=4)         :: my_number
     integer(kind=4)         :: my_timestep
     real(kind=4)            :: my_form_aexp ! expansion factor when halo got half its mass
     integer(kind=4)         :: level, hosthalo, hostsub,nbsub,nextsub ! data for structure tree
     real(kind=4)            :: m

     real(kind=4)            :: r
     real(kind=4)            :: spin
     real(kind=4)            :: ek,ep,et
  end type halo

 integer(kind=4)                 :: numero_step

end module tree_defs

Program main

  use ramses_util
  use tree_defs
  implicit none

  real(kind=8), parameter :: Msol = 1.9891d+33
  real(kind=8), parameter :: Mparsec = 3.08d+24

  type(halo),allocatable,target   :: halo_list(:) ! list of all halos 
  type(halo)   :: h
  integer :: nb_of_halos, nb_of_subhalos,nbodies,nbtot,jj,ji,jf,ind,k,npart,nx
  integer :: j,unitfile,i,lounode,lounei
  character(128)   :: treefile,filenode,filenei
  real(kind=4), dimension(3) :: pos
  integer(kind=4)             :: nparts_inhalo
  character(80)               :: filename,file
  integer(kind=4),allocatable :: indices(:)
  integer(kind=4),allocatable :: strt(:,:)
  integer(kind=4),allocatable :: strtp(:),strti(:),struct(:)
  real(kind=4)                :: massp,valmax,omega_t,age_univ,aexp0
  real(kind=8)                :: scale_l,scale_d,scale_t, &
                                 aexp,omega_m,omega_l,omega_k,h0,lbox,t,boxlen
  integer(kind=4),allocatable     :: liste_parts(:)
  character(len=128) :: arg,repository
  logical::ok

  call getarg(2,arg)
  repository = trim(arg)

  call read_info(arg,boxlen,scale_l,scale_d,scale_t,t,aexp,omega_m,omega_l,omega_k)

  print *,'to find the size of the box (Mpc) '
  print *,scale_l,scale_d,scale_t,aexp
  print *,'boxlen ',scale_l/aexp/Mparsec,scale_l/Mparsec,scale_l/aexp/Mparsec*0.7

  lbox=scale_l/Mparsec

!  filenode='halomaker.nod'
  call getarg(3,arg)
  filenode=trim(arg)

  call getarg(4,arg)  
  filenei=trim(arg)

  lounode=20
  open(unit=lounode,status='unknown',form='formatted',file=filenode)

  unitfile = 10
  call getarg(5,arg)
  filename = trim(arg)
  print *,'filename ',filename
  !    filename='tree_bricks044'
  !! Leo: check if tree_brick file exists
  inquire(file=filename, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(filename)//' not found.'
     stop
  endif
  open(unit=unitfile,file=filename,form='unformatted',status='old')
  read(unitfile) nbodies
  print *,'nbodies ',nbodies

  read(unitfile) massp
  print *,'massp ',massp

  read(unitfile) aexp0
  read(unitfile) omega_t
  read(unitfile) age_univ
  read(unitfile) nb_of_halos, nb_of_subhalos
  
  print *,'aexp, omega_t, age_univ ',aexp0, omega_t, age_univ
  print *,'nb_of_halos ',nb_of_halos, nb_of_subhalos

  nbtot=nb_of_halos+nb_of_subhalos
  print *,'nbtot ',nbtot

  allocate(halo_list(nbtot))
  allocate(strtp(nbodies))
  allocate(strti(nbodies))
  allocate(struct(nbodies))

  struct=0
  valmax=0.
  jj=0
  ji=1
  do i=1,nbtot
     read(unitfile) nparts_inhalo
     jj=jj+nparts_inhalo
     allocate(indices(nparts_inhalo))
     read(unitfile) indices

     do k=1,nparts_inhalo
        ind=indices(k)
        struct(ind)=i
     enddo

     call read_halo(halo_list(i),unitfile)
     pos(1)=halo_list(i)%p%x
     pos(2)=halo_list(i)%p%y
     pos(3)=halo_list(i)%p%z

     if ( (i<10) .or. (halo_list(i)%my_number == 688) ) print *, pos(1:3)

     pos(1)=pos(1)/lbox+.5
     pos(2)=pos(2)/lbox+.5
     pos(3)=pos(3)/lbox+.5

     if ( (pos(1)>1.) .or. (pos(2)>1.) .or. (pos(3)>1.) ) print *,'stop ',pos(1:3)
     if ( (pos(1)>1.) .or. (pos(2)>1.) .or. (pos(3)>1.) ) stop

     halo_list(i)%r=halo_list(i)%r/lbox
     halo_list(i)%datas%rvir=halo_list(i)%datas%rvir/lbox

     valmax=max(valmax,halo_list(i)%datas%rvir)

     if ( (i<10) .or. (halo_list(i)%my_number == 688) ) then
        print *,i,nparts_inhalo
        print *,pos(1:3)
        print *,halo_list(i)%my_number,halo_list(i)%level,halo_list(i)%nbsub
        print *,halo_list(i)%hosthalo,halo_list(i)%hostsub
        print *,halo_list(i)%m*1.e11,halo_list(i)%p%x,halo_list(i)%p%y,halo_list(i)%p%z
        print *,halo_list(i)%r,halo_list(i)%datas%rvir,halo_list(i)%datas%mvir*1.e11
        print *
     endif

     halo_list(i)%m=halo_list(i)%m*1.e11
     halo_list(i)%datas%mvir=halo_list(i)%datas%mvir*1.e11
     
     write(lounode,'(I,I,I,I,I,I,E,E,E,E,I,E,E,E,E)') &
          &                         halo_list(i)%my_number, &
          &                         halo_list(i)%level, &
          &                         halo_list(i)%hosthalo, &
          &                         halo_list(i)%hostsub, &
          &                         halo_list(i)%nbsub, &
          &                         halo_list(i)%nextsub, &
          &                         halo_list(i)%et, &
          &                         halo_list(i)%datas%mvir, &
          &                         halo_list(i)%datas%rvir, &
          &                         halo_list(i)%r, &
          &                         nparts_inhalo, &
          &                         halo_list(i)%m, &
          &                         pos(1:3)
     
     jf=ji+nparts_inhalo-1
     strti(ji:jf)=halo_list(i)%my_number
     strtp(ji:jf)=indices(:nparts_inhalo)
     ji=jf+1

     deallocate(indices)
  enddo
  print *,'END loop'
  close(unitfile)
  close(lounode)

  print *,'valmax ',valmax
  print *,'jj ',jj,jf

  lounei=22
  open(unit=lounei,status='unknown',form='unformatted',file=filenei)   
  write(lounei) nbodies
  write(lounei) struct(:nbodies)
  close(lounei)
  
  deallocate(halo_list)
  deallocate(strti)
  deallocate(strtp)

end program main

!***********************************************************************
subroutine read_halo(h,unitfile)

  ! reads halo properties computed by part_to_halo.F

  use tree_defs
  implicit none

  integer(kind=4) :: unitfile
  type (halo)     :: h

  read(unitfile) h%my_number
  read(unitfile) h%my_timestep
  ! change the time step number to match the number of branching (time resolution) 
  ! you decided your tree is going to have 
  !    h%my_timestep = numero_step
  read(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
  read(unitfile) h%m
  read(unitfile) h%p%x,h%p%y,h%p%z
  read(unitfile) h%v%x,h%v%y,h%v%z
  read(unitfile) h%L%x,h%L%y,h%L%z 
  read(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
  read(unitfile) h%ek,h%ep,h%et
  read(unitfile) h%spin
  read(unitfile) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
  read(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c

  return

end subroutine read_halo
!***********************************************************************
