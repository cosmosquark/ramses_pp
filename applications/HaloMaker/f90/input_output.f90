module input_output

  !! JB: This version only deals with outputs form RAMSES (v 3.07) or Gadget

  use subbox

  public

contains

!///////////////////////////////////////////////////////////////////////
!***********************************************************************
  subroutine read_data

    ! This routine read the output of N-body simulations (particles positions and speeds, 
    ! cosmological and technical parameters)

    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ! WARNING: this routine just reads the data and converts positions       !
    !          and velocities from CODE units to these units                 !
    !          -- positions are between -0.5 and 0.5                         !
    !          -- velocities are in km/s                                     !
    !             in units of Hubble velocity accross the box for SIMPLE (SN)!
    !          -- total box mass is 1.0                                      !
    !             for simulation with hydro (only with -DRENORM) flag        !
    !          -- initial (beg of simulation) expansion factor is ai=1.0     ! 
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    implicit none

    character(len=200) :: name_of_file

    write(errunit,*)
    write(errunit,*) '> In read_data: timestep  ---> ',numero_step

    if (numero_step == 1) then  
       write(errunit,*) '> data_dir: ''',trim(data_dir),''''
       ! contains the number of snapshots to analyze and their names, type and number (see below)
       write(name_of_file,'(a,a)') trim(data_dir),'inputfiles_HaloMaker.dat'
       open(unit=12,file=name_of_file,status='old')
       !call open_exist(12,'inputfiles_HaloMaker.dat')
    endif

    ! then read name of snapshot, its type (pm, p3m, SN, Nzo, Gd), num of procs used and number of snapshot
    read(12,*) name_of_file,type,nbPes,numstep
    write(file_num,'(i3.3)') numstep
    ! jeje 
    !    write(name_of_file,'(a,a)') trim(data_dir),trim(name_of_file)
    ! end jeje 

    if (numero_step == nsteps) close(12)

    if (type == 'Ra3') then

       if (subboxFoF) then 
#ifdef STARS 
          call read_ramses_new_stars_subbox(name_of_file) 
#else
          call read_ramses_new_subbox(name_of_file)
#endif
       else
#ifdef STARS
          call read_ramses_new_stars(name_of_file)
#else
          call read_ramses_new(name_of_file)
#endif
          nbodies_full_sim = nbodies
       end if

       ! Computation of omega_t = omega_matter(t)
       !
       !                            omega_f*(1+z)^3
       ! omega(z)   = ------------------------------------------------------------------
       !              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f
       !
       !
       !                              omega_lambda_0
       ! omega_L(z) = ----------------------------------------------------------------
       !              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f

       omega_t  = omega_f*(af/aexp)**3
       omega_t  = omega_t/(omega_t+(1.-omega_f-omega_lambda_f)*(af/aexp)**2+omega_lambda_f)

    else if (type == 'Gd') then
       
       if (subboxFoF) then
          call read_gadget_subbox(name_of_file)
       else
          call read_gadget(name_of_file)
          ! jeje 
          nbodies_full_sim = nbodies
          ! end jeje
       end if
       
       ! Computation of omega_t = omega_matter(t)
       !
       !                            omega_f*(1+z)^3
       ! omega(z)   = ------------------------------------------------------------------
       !              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f
       !
       !
       !                              omega_lambda_0
       ! omega_L(z) = ----------------------------------------------------------------
       !              omega_f*(1+z)^3+(1-omega_f-omega_lambda_f)*(1+z)^2+omega_lambda_f

       omega_t  = omega_f*(af/aexp)**3
       omega_t  = omega_t/(omega_t+(1.-omega_f-omega_lambda_f)*(af/aexp)**2+omega_lambda_f)
       
    else
       write(errunit,*) '> Don''t know the snapshot format: ''',type,''''
       stop
    endif
    
    write(errunit,*) '> min max position (in box units)   :',minval(pos),maxval(pos)
    write(errunit,*) '> min max velocities (in km/s)      :',minval(vel),maxval(vel)
    write(errunit,*) '> Reading done.'
    write(errunit,*) '> aexp = ',aexp

    return

  end subroutine read_data

!***********************************************************************
  subroutine open_exist(unit,name)

    implicit none

    integer(kind=4)                                  :: unit,lnblnk
    character(len=*)                                 :: name
    character(len=len_trim(name)+len_trim(data_dir)) :: file

    write(file,'(a,a)') data_dir(1:lnblnk(data_dir)),name(1:lnblnk(name))
    open(unit=unit,file=file,status='unknown')

    return

  end subroutine open_exist

!***********************************************************************
  subroutine open_append(unit,name)

    implicit none

    integer(kind=4)                                  :: unit,lnblnk
    character(len=*)                                 :: name    
    character(len=len_trim(name)+len_trim(data_dir)) :: file

    write(file,'(a,a)') data_dir(1:lnblnk(data_dir)),name(1:lnblnk(name))
    if (numero_step == 1) then 
       open(unit=unit,file=file,status='unknown')
    else
       open(unit=unit,file=file,status='old',position='append')
    end if

    return

  end subroutine open_append

!***********************************************************************
  subroutine read_gadget(name_file)


    ! Useful info for reading Gadget2
    ! Particles in gadget are sorted into 6 type
    ! The number of particles in each type is found in io_header%Vnpart_tot(1:6) and their masses in io_header%massarr(1:6)
    ! > When output is splitted the number of particles per file is io_header%Vnpart(1:6)
    ! > When output is splitted the number of particles is io_header%Vnpart(1:6) = io_header%Vnpart_tot(1:6) 

    ! > In most cases
    ! particles type 1 Gas
    ! particles type 2 Dark Matter
    ! particles type 5 Star
    ! nbodies = io_header%Vnpart_tot(2) 
    ! nbaryon = io_header%Vnpart_tot(1) + io_header%Vnpart_tot(5) 
    ! id of DM particles are between nbaryon - 1 and nbaryon + nbodies - 1

    ! > For simulation -DBIG_RUN option
    ! io_header%massarr(2) > 0., Gd_mass is not written in the file
    ! > For resimulation 
    ! io_header%massarr(2) = 0., Gd_mass is written in the file
   
    ! > If resimulation dark matter only ( Gadget 1 format )
    ! particles type 2 high resolution Dark Matter particles
    ! particles type 3 low  resolution Dark Matter particles
    ! particles type 4 low  resolution Dark Matter particles
    ! particles type 5 low  resolution Dark Matter particles
    ! nbodies = io_header%Vnpart_tot(2) + io_header%Vnpart_tot(3) + io_header%Vnpart_tot(4) + io_header%Vnpart_tot(5) 
    ! nhr     = io_header%Vnpart_tot(2)
    ! nlr     = io_header%Vnpart_tot(2) + io_header%Vnpart_tot(3) + io_header%Vnpart_tot(4) + io_header%Vnpart_tot(5) 
    ! nbaryon = io_header%Vnpart_tot(2) + io_header%Vnpart_tot(1) + io_header%Vnpart_tot(6) 
    ! All values of io_header%massarr are 0, all other values of io_header%Vnpart_tot should be 0
    ! id of DM particles are between 1 and nbodies 


    implicit none

    integer(kind=4)                      :: nbig,stop_read,nf,istat
    character(len=*)                     :: name_file
    character(len=3)                     :: fnumber
    character(len=len_trim(name_file)+4) :: Vfilename
    type header
       integer(kind=4) :: Vnpart(6),Vnpart_tot(6)
       real(kind=8)    :: massarr(6),npartTotalHighWord(6)
       real(kind=8)    :: time,redshift,Boxsize,Omega0,OmegaLambda,HubbleParam
       integer(kind=4) :: flag_sfr,flag_cooling,flag_feedback,flag_stellarage,flag_metals,flag_entropy_instead_u
       integer(kind=4) :: num_files
!       character*60    :: fill
    end type header
    type(header)       :: io_header
    integer(kind=4)                      :: VN,nbaryon,id_part,idp,itype
    integer(kind=4),allocatable          :: Vid(:)
    real(kind=4),allocatable             :: Gd_pos(:,:),Gd_vel(:,:)    
#ifndef BIG_RUN
    real(kind=4),allocatable             :: Gd_mass(:)
    integer(kind=4)                      :: VN_lr,VN_total,nlr
    real(kind=4)                         :: HRsoftlen, LRsoftlen    
#endif
    real(kind=4)                         :: Lbox_pt_local
    ! softening lengths (comoving, and maximum physical in kpc/h) used in the simulation:
    real(kind=4),parameter               :: HRepscom = 10.0, HRepsmaxph = 5.0
    real(kind=4),parameter               :: LRepscom = 10.0, LRepsmaxph = 5.0
    integer(kind=4)                      :: i,ierr
    real(kind=4)                         :: Omega0,OmegaBaryon,mass_fac,mass_box 
    character*200                        :: paramname,value,line,fileparam
    logical(kind=4)                      :: gadget1
#ifdef Test_Gd
    integer(kind=4)                      :: typeoffset,ngood,nbad,nrest
    integer(kind=4),allocatable          :: typeindex(:)
#ifndef BIG_RUN
    integer(kind=4)                      :: nvlr
    real(kind=4)                         :: massmin,massmax
    real(kind=4), allocatable            :: massloc(:)
#endif
#endif

    ! first build file name as path/snapshot_xxx.n, where:
    ! "path/snapshot_xxx" are in "name_file" variable, and 
    ! "n" is built by: //fnumber(verify(fnumber,' '):3)
    stop_read = 0
    nf        = -1
    nbig      = 1000
    nbodies   = 0
    nbaryon   = 0
    gadget1   = .false.

!********************************************************
    read_files: do while (nf <= nbig .and. stop_read == 0)
       istat = 1
       ! serch the file to read
       if(nf.lt.0) then
          ! try to open the single file if it exist
          Vfilename = trim(name_file)
          open(1,file=Vfilename,status='old',form='unformatted',iostat=istat)
          if(istat == 0) then
             stop_read = 1
             write(errunit,*) '> read_gadget... only one file to read'
          else
             ! .. if something went wrong (iostat/=0): 
             write(errunit,*) '> read_gadget... several files to read'
          end if
       end if
       if(istat /= 0) then
          nf        = nf+1
          write(fnumber,'(i3)') nf  !--> transform the integer nf in a character
          Vfilename = trim(name_file)// '.'//fnumber(verify(fnumber,' '):3)
          ! .. try to open file # (nf-1)
          ! possibly only several snap per step exist, thus .fnumber exists 
          ! (the file name is snapshot_xxx.fnumber instead of snapshot_xxx), so 
          ! try to read again, wit .fnumber: starting at zero if .1 isn't found either we suppose no file exist
          open(1,file=Vfilename,status='old',form='unformatted',iostat=istat)
          ! .. if something went wrong again ==> stop here:
          if (istat /= 0.and.nf.le.1) then
             write(errunit,*) ' > ERROR: read_gadget: input file not there',trim(Vfilename)
             stop
          end if
          if(istat /=0) stop_read = 1
       endif
       if (istat /= 0) exit read_files

       write(errunit,*) '> Read file: ',trim(Vfilename)

       ! .. if everything ok, continue:
       read (1,iostat=istat)                  &
            io_header%Vnpart,                 &  
            io_header%massarr,                &  
            io_header%time,                   &
            io_header%redshift,               &
            io_header%flag_sfr,               &
            io_header%flag_feedback,          &
            io_header%Vnpart_tot,             &
            io_header%flag_cooling,           &
            io_header%num_files,              &
            io_header%Boxsize,                &
            io_header%Omega0,                 &
            io_header%OmegaLambda,            &
            io_header%HubbleParam,            &
            io_header%flag_stellarage,        &
            io_header%flag_metals,            &
            io_header%npartTotalhighWord,     &
            io_header%flag_entropy_instead_u
       ! .. if file exists but is empty, stop reading
       if (istat /= 0) stop ' > ERROR: read_gadget... could not read header.'

       if(nbodies.le.0) then
          ! if nbodies is not allocated we're reading the first file
          
          write(errunit,*)
          write(errunit,*) '> Useful data from header:'
          write(errunit,*) '> io_header%time (aexp/af)          :',real(io_header%time,4)
          write(errunit,*) '> io_header%redshift                :',real(io_header%redshift,4)
          write(errunit,*) '> io_header%Boxsize                 :',real(io_header%Boxsize,4)
          write(errunit,*) '> io_header%Omega0                  :',real(io_header%Omega0,4)
          write(errunit,*) '> io_header%OmegaLambda             :',real(io_header%OmegaLambda,4)
          write(errunit,*) '> io_header%HubbleParam             :',real(io_header%HubbleParam,4)
          write(errunit,*)

#ifdef Test_Gd
          do itype = 1,6
             write(errunit,'(1x,a,i1,1x,a1,1x,i10)') '> Tot nb of particles type ',itype,':',io_header%Vnpart_tot(itype)
             write(errunit,'(1x,a,i1,1x,a1,1x,e14.6)') '> mass for particles type  ',itype,':',real(io_header%massarr(itype))
          end do
          write(errunit,*)
#endif
          nbodies = io_header%Vnpart_tot(2) 
          nbaryon = io_header%Vnpart_tot(1) + io_header%Vnpart_tot(5)
          if(io_header%Vnpart_tot(3).gt.0.or.io_header%Vnpart_tot(4).gt.0.or.io_header%Vnpart_tot(6).gt.0) then
             gadget1 = .true.
#ifdef BIG_RUN       
             write(errunit,*)
             write(errunit,*) '> WARNING: the code was compiled with -DBIG_RUN option'
             write(errunit,*) '> Some particles have not been accounted for.'
             do itype = 1,6
                write(errunit,'(1x,a,i1,1x,a1,1x,i10)') '> Tot nb of particles type ',itype,':',io_header%Vnpart_tot(itype)
                write(errunit,'(1x,a,i1,1x,a1,1x,e14.6)') '> mass for particles type  ',itype,':',real(io_header%massarr(itype))
             end do
             write(errunit,*) '> Make sure this is correct'
             write(errunit,*)
             if(io_header%Vnpart(2).gt.0.and.io_header%massarr(2).eq.0.) then
                write(errunit,*)
                write(errunit,*) '> ERROR: read_gadget... mass is 0. for DM particles (type2)'
                write(errunit,*) '> Code was compiled with -DBIG_RUN option'
                write(errunit,*) '> This shouldn''t be the case'
                stop
             end if
#else
             do itype = 1,6
                write(errunit,'(1x,a,i1,1x,a1,1x,i10)') '> Tot nb of particles type ',itype,':',io_header%Vnpart_tot(itype)
                write(errunit,'(1x,a,i1,1x,a1,1x,e14.6)') '> mass for particles type  ',itype,':',real(io_header%massarr(itype))
             end do
             write(errunit,*)
             write(errunit,*) '> WARNING: interpreting type 2 particles as high res DM particles '
             nhr = io_header%Vnpart_tot(2)
             if(io_header%Vnpart(1).gt.0) then
                write(errunit,*) '> WARNING: interpreting type 3, type 4 particles as low res DM particles'
                nlr = io_header%Vnpart_tot(3) + io_header%Vnpart_tot(4)
                write(errunit,*) '> WARNING: interpreting type 1, type 5 particles as baryonnic matter'
                write(errunit,*) 
             else
                write(errunit,*) '> WARNING: interpreting type 3, type 4,type 5 particles as low res DM particles'
                nlr = io_header%Vnpart_tot(3) + io_header%Vnpart_tot(4) + io_header%Vnpart_tot(5)
             endif
             if(io_header%Vnpart_tot(6).gt.0) write(errunit,*) '> WARNING: don''t know how to interpret type 6 particles'
             write(errunit,*) '> Rerun after compiling with -DTest_Gd option for testing'
             write(errunit,*) '> Nb of high resolution particles   :',nhr 
             write(errunit,*) '> Nb of low resolution particles    :',nlr
             nbodies = nhr + nlr
             nbaryon = io_header%Vnpart_tot(1) + io_header%Vnpart_tot(6)
#endif
          endif
          write(errunit,*) '> Number of DM particles            :',nbodies
          write(errunit,*) '> Number of baryons                 :',nbaryon
            
          allocate(pos(nbodies,3))
          allocate(vel(nbodies,3))
          ! jeje 
          allocate(ParticleID(nbodies))
          particleID = -1
          ! end jeje 
          pos     = 0.0
          vel     = 0.0
#ifndef BIG_RUN
          allocate(mass(nbodies), epsvect(nbodies))
          mass    = 0.0
          epsvect = 0.0
#endif
       else
          gadget1 = .false.
       end if

#ifdef Test_Gd
       if(nf.ge.0) then
          ! For now we don't know to what kind of particles type 3,4 and 6 do correspond to
          write(errunit,*)
          do itype = 1,6
             write(errunit,'(1x,a,i1,1x,a1,i10)') '> Nb of particles type ',itype,':',io_header%Vnpart(itype)
          end do
          write(errunit,*)
       end if
#endif

       VN    = sum(io_header%Vnpart)       
       write(errunit,*) '> Total nb of particles in file     :',sum(io_header%Vnpart)

       allocate(Gd_pos(3,VN),Gd_vel(3,VN),Vid(VN))
       read(1) Gd_pos
       read(1) Gd_vel
       read(1) Vid
#ifndef BIG_RUN
       allocate(Gd_mass(VN))
       if(gadget1.and.io_header%massarr(2).gt.0..and.io_header%Vnpart(1).eq.0) then
           read(1,iostat=istat) Gd_mass(io_header%Vnpart(2)+1:) !old version Hr masses were not written it seems.
       else
          read(1,iostat=istat) Gd_mass
       end if
       if(istat /= 0) stop ' > ERROR: read_gadget... masses are not written in this file.'
#endif
       close(1)

#ifdef Test_Gd
       ! We thing that DM part id are written with an offset corresponding to nbaryons
       write(errunit,'(1x,a,6(i10,1x))')   '> io_header%Vnpart :',io_header%Vnpart
       write(errunit,'(1x,a,6(E10.3,1x))') '> io_header%massarr:',io_header%massarr
       write(errunit,*)
       typeoffset = 0
       do itype = 1, 6
          write(errunit,'(1x,a,i1,a,i10)') '> type: ',itype,' io_header%Vnpart ',io_header%Vnpart(itype)
          if(io_header%Vnpart(itype).gt.0) then
             write(errunit,*) '> between:', typeoffset +1,'and', typeoffset+io_header%Vnpart(itype)
             ngood = 0
             nbad  = 0
             nrest = 0
             if(nbaryon .gt. 0) then
                do id_part = typeoffset+1, typeoffset+io_header%Vnpart(itype)
                   
                   if(Vid(id_part).ge.0.and.Vid(id_part).lt.nbaryon) then
                      ngood = ngood + 1
                   else if(Vid(id_part).lt.nbaryon + nbodies) then
                      nbad  = nbad + 1
                   else
                      nrest = nrest + 1
                   end if
                end do
                write(errunit,*) '> id between 0 and nbaryon-1               :',ngood
                write(errunit,*) '> id between nbaryon and nbaryon+nbodies-1 :',nbad
                write(errunit,*) '> id below 0 or above nbaryon+nbodies-1    :',nrest   
             else
                do id_part = typeoffset+1, typeoffset+io_header%Vnpart(itype)
                   if(Vid(id_part).gt.0.and.Vid(id_part).le.nbodies) then
                      ngood = ngood + 1
                   else if(Vid(id_part).eq.0) then
                      nbad  = nbad + 1
                   else
                      nrest = nrest + 1
                   end if
                end do
                write(errunit,*) '> id between 1 and nbodies                 :',ngood
                write(errunit,*) '> id between eq 0                          :',nbad
                write(errunit,*) '> id below 0 or above nbodies              :',nrest                   
             end if
             allocate(typeindex(io_header%Vnpart(itype)))
             typeindex(1:io_header%Vnpart(itype)) = Vid(typeoffset+1:typeoffset + io_header%Vnpart(itype))
             write(errunit,*) '> main,max Vid',minval(typeindex),maxval(typeindex)
             deallocate(typeindex)
#ifndef BIG_RUN
             allocate(massloc(io_header%Vnpart(itype)))
             massloc(1:io_header%Vnpart(itype)) = Gd_mass(typeoffset+1:typeoffset + io_header%Vnpart(itype))
             write(errunit,*) '> main,max Gd_mass',minval(massloc),maxval(massloc)
             deallocate(massloc)
#endif
          end if
          typeoffset = typeoffset + io_header%Vnpart(itype)
          write(errunit,*)
       end do
#endif
       
       ! if SPH simulation, you have gas particles (number in io_header%Vnpart(1))
       ! which can turn into star particles (in io_header%Vnpart(5)) so we have to skip
       ! these to read the DM particle data only
       
       do id_part = io_header%Vnpart(1)+1,io_header%Vnpart(1)+io_header%Vnpart(2) ! loop only on DM particles 
          idp   = Vid(id_part) - nbaryon
          if (idp == 0) idp = nbodies ! offset from C to FORTRAN numbering
          if (idp < 1 .or. idp > nbodies ) then
             write(errunit,*) '> idp,nbodies',idp,nbodies
             stop ' > ERROR: read_gadget... bad particle id'
          end if
          pos(idp,:) = Gd_pos(:,id_part)     
          vel(idp,:) = Gd_vel(:,id_part)
#ifndef BIG_RUN
          mass(idp)  = Gd_mass(id_part)
#endif
          ! jeje
          particleID(idp) = VID(id_part)
       end do
       
#ifndef BIG_RUN
       if(gadget1) then
          !loop for low rest of particles
          do id_part = io_header%Vnpart(1)+io_header%Vnpart(2) +1,&
               io_header%Vnpart(1)+io_header%Vnpart(2) +io_header%Vnpart(3)+io_header%Vnpart(4)
             idp   = Vid(id_part) - nbaryon
             if (idp == 0) idp = nbodies ! offset from C to FORTRAN numbering
             if (idp < 1 .or. idp > nbodies ) then
                write(errunit,*) '> idp,nbodies',idp,nbodies
                stop ' > ERROR: read_gadget... bad particle id'
             end if
             pos(idp,:) = Gd_pos(:,id_part)     
             vel(idp,:) = Gd_vel(:,id_part)
             mass(idp)  = Gd_mass(id_part)
          end do
          if(io_header%Vnpart_tot(1).eq.0) then
             do id_part =   io_header%Vnpart(1)+io_header%Vnpart(2) +io_header%Vnpart(3)+io_header%Vnpart(4)+1,&
                  io_header%Vnpart(1)+io_header%Vnpart(2) +io_header%Vnpart(3)+io_header%Vnpart(4)+io_header%Vnpart(5)
                idp   = Vid(id_part) - nbaryon
                if (idp == 0) idp = nbodies ! offset from C to FORTRAN numbering
                if (idp < 1 .or. idp > nbodies ) then
                   write(errunit,*) '> idp,nbodies',idp,nbodies
                   stop ' > ERROR: read_gadget... bad particle id'
                end if
                
                pos(idp,:) = Gd_pos(:,id_part)     
                vel(idp,:) = Gd_vel(:,id_part)
                mass(idp)  = Gd_mass(id_part)
             end do
          end if
       end if
#endif
       
#ifndef BIG_RUN
       deallocate(Gd_mass)
#endif
       deallocate(Gd_pos, Gd_vel, Vid)
       
    enddo read_files
!********************************************************
    
    aexp  = io_header%time*af  ! because SN format assumes a(t_in)=1
!    write(errunit,*) '  SN_aexp, z                                   : ',aexp,-1.+af/aexp 
#ifdef BIG_RUN
    ! have to recalculate mboxp because Gadget particles are in physical units 
    mass_box = io_header%massarr(2)*10./H_f*real(nbodies,4) ! in 10^11 M_sun
    massp    = 1./real(nbodies,4) !in units of mboxp !io_header%massarr(2)*10./H_f/mass_box  
#else
    ! count low res and high res DM particles
    massp = minval(mass)
    nhr     = 0
    nlr     = 0
    do idp = 1, nbodies
       if(mass(idp).eq.massp) then
          nhr = nhr + 1
       else if(mass(idp).gt.massp) then
          nlr = nlr + 1
       end if
    end do
    write(errunit,*) '> Nb of high resolution particles   :',nhr 
    write(errunit,*) '> Nb of low resolution particles    :',nlr
    
    ! .. check that reading went ok for resimulation
    VN_total = nhr+nlr
    Vn_lr    = nlr
    if (VN_total /= nbodies) then
       stop ' > read_gadget... ERROR: nlr + nhr /= nbodies !' 
    endif
    mass_box = real(sum(real(mass,8)),4)*10./H_f ! in 10^11 M_sun
    massp    = massp*10./H_f/mboxp  !in units of mboxp
#endif
    
    ! Change units of pos, vel, masses: 
    ! .. Transform comoving positions in h^-1 Kpc to physical positions in Mpc, and 
    !    normalise to physical box length in Mpc to have them in [-0.5,0.5]:
    !    Gd_pos * h^-1 * 10^-3 *io_header%time/Lbox(t), where Lbox(t)=Lboxp*aexp
    !    pos  = pos*io_header%time/(10.*H_f*Lboxp*aexp) - 0.5
    !    We can as well normalize with io_header%Boxsize = 10.*H_f*Lboxp
    ! .. Transform Gadget velocities, in SIMPLE_peculiar_velocies in units
    !    of Hubble parameter accross the box  
    !  where:
    !     Gd_vel= sqrt(io_header%time)*x_dot in km/s [io_header%time == Gadget expansion param]
    !  thus: 
    !     pec_vel=[sqrt(io_header%time)*Gd_vel/H(t)/L(t)]

    pos           = pos*io_header%time/(io_header%Boxsize*aexp) - 0.5
    Lbox_pt_local = Lboxp*(aexp/ai)
    vel           = vel *sqrt(io_header%time)
#ifndef BIG_RUN
    ! Physical softening length at the current time, used in the original simulation
    HRsoftlen  = min(real(io_header%time,4)*HRepscom,HRepsmaxph)
    LRsoftlen  = min(real(io_header%time,4)*LRepscom,LRepsmaxph)
    ! get epsvect = softening length in units of "fraction of Lbox(t)", as needed by subroutine softgrav:
    epsvect(1:nhr)  = HRsoftlen/(10.*H_f*Lbox_pt_local)
    epsvect(nhr+1:) = LRsoftlen/(10.*H_f*Lbox_pt_local)
    !write(errunit,*) ' '
    !write(errunit,*) 'softening length of HR parts in units of Lbox(t)=', epsvect(1:1)
    !write(errunit,*) 'softening length of LR parts in units of Lbox(t)=', epsvect(nhr+1:nhr+1)
    ! .. Transform masses from 10^10 Msol/h --> to units of mboxp
    mass          = mass * 10.0 / H_f / mass_box !in units of box total mass
    ! .. min and max mass:
    !write(errunit,*) 
    minlrmrat = minval(mass(nhr+1:))/massp
    !write(errunit,*) '  min-max MASS of LR particles [in units of mhr]  : ', &
    !                    minlrmrat,maxval(mass(nhr+1:))/massp
#endif  

    if (nbaryon > 0) then 
       write(errunit,*)
       write(errunit,*) '> Checking the "parameters-usedvalues" file: '
       i         = index(name_file,'/snapshot_')
       fileparam = trim(name_file(1:i))//'parameters-usedvalues'
       open(unit=5,status='old',form='formatted',file=fileparam,iostat=ierr)
       if(ierr.ne.0) then
          write(errunit,*) '> Please copy or link the "parameters-usedvalues" file in directory:'
          write(errunit,*) fileparam(1:i)
          stop
       else
          Omega0       = -1.
          OmegaBaryon  = -1.0
          do
             read(5,'(a)',end=51) line
             i = scan(line,'     ')
             if(i.eq.0.or.line(1:1) .eq. '#') cycle
             paramname = trim(adjustl(line(:i-1)))
             value     = trim(adjustl(line(i+1:)))
             select case(trim(paramname))
             case('Omega0')
                read(value,*) Omega0
             case('OmegaBaryon')
                read(value,*) OmegaBaryon
             end select
          end do
51        close(5)
          write(errunit,*) '> paremeters I neaded form the "parameter-usedvalues" file'
          write(errunit,*) '> Omega0       : ',Omega0
          write(errunit,*) '> OmegaBaryon  : ',OmegaBaryon
          if (Omega0 .le. 0.) stop ' > read_gadget... Coudn''t read Omega0 value'
          if (OmegaBaryon .le. 0.) stop ' > read_gadget... Coudn''t read OmegaBaryon value'
          write(errunit,*)
       end if
       mass_fac = 1./(1.0-OmegaBaryon/Omega0)
       mass_box = mass_box * mass_fac
#ifndef RENORM
       massp = massp/mass_fac
#ifndef BIG_RUN
       mass  = mass/mass_fac
#endif
#endif
    endif

    write(errunit,*) '> In Gadget'
    write(errunit,*) '> box mass (M_sun)            :',mass_box*1.e11
    write(errunit,*) '> particle mass (in M_sun)    :',massp*mass_box*1.e11         
#ifndef BIG_RUN
    write(errunit,*) '> HR particle mass (in M_sun) :', minval(mass)*mass_box*1e11
#endif
#ifdef RENORM
    if(nbaryon > 0) then
       write(errunit,*) '> after renormalisation'
       write(errunit,*) '> box mass (M_sun)         : ',mass_box*1e11
       write(errunit,*) '> particle mass (M_sun)    :', massp*mass_box*1e11
    end if
#endif
    write(errunit,*)
    write(errunit,*) '> In HaloMaker, after rescaling'
    write(errunit,*) '> box mass (M_sun)            :',mboxp*1.e11
    write(errunit,*) '> particle mass (in M_sun)    :',massp*mboxp*1.e11         
#ifndef BIG_RUN
    write(errunit,*) '> HR particle mass (in M_sun) :', minval(mass)*mass_box*1e11
#endif
#ifdef RENORM
    if(nbaryon > 0) then
       write(errunit,*) '> after renormalisation'
       write(errunit,*) '> box mass (M_sun)         : ',mboxp*1e11
       write(errunit,*) '> particle mass (M_sun)    :', massp*mboxp*1e11
    end if
#endif
 

    return

  end subroutine read_gadget

!****************************************************************************************************
  subroutine read_ramses_new_stars(repository)
    
    ! This routine reads STELLAR particles dumped in the RAMSES 3.0 new i/o format.  

    ! repository is directory containing output files
    ! e.g. /horizon1/teyssier/ramses_simu/boxlen100_n256/output_00001/

    implicit none

    character(len=*)            :: repository
    integer(kind=4)             :: ndim,npart,idim,icpu,ipos,ncpu,i,ipar
    integer(kind=4)             :: ncpu2,npart2,ndim2,idum,nout,nsink,nstar
    integer(kind=4)             :: nx,ny,nz,nlevelmax,ngridmax,nstep_coarse
    integer(kind=4),allocatable :: levelp(:)
    integer(kind=4),allocatable :: idp(:)
    real(kind=8)                :: boxlen,tco,aexp_ram,hexp
    real(kind=8)                :: omega_m,omega_l,omega_k,omega_b
    real(kind=8)                :: scale_l,scale_d,scale_t,dummy
    real(kind=8)                :: massres
    real(kind=8),allocatable    :: dumout(:),tmpp(:,:),tmpv(:,:),tmpm(:),tmpt(:)
    character(len=200)          :: line,name,value,nomfich
    character(len=5)            :: nchar,ncharcpu
    logical                     :: ok
    integer(kind=4) :: cnt
    real(kind=8),parameter :: gramm_to_1011Msun = 1.d0 / 1.9891d+33 / 1.d11

    ! read cosmological params in header of amr file
    ipos    = index(repository,'output_')
    nchar   = repository(ipos+7:ipos+13)
    nomfich = trim(repository)//'amr_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok)
    if (.not. ok) then
       write(errunit,*)'File '//trim(nomfich)//' not found'
       stop
    else
       open(unit=10,file=nomfich,status='old',form='unformatted')
       read(10) ncpu
       read(10) ndim
       read(10) nx,ny,nz
       read(10) nlevelmax
       read(10) ngridmax
       read(10) idum
       read(10) idum
       read(10) boxlen
       read(10) nout,idum,idum
       allocate(dumout(nout))
       read(10) dumout
       read(10) dumout
       deallocate(dumout)
       read(10) tco
       allocate(dumout(nlevelmax))
       read(10) dumout
       read(10) dumout
       deallocate(dumout)
       read(10) idum,nstep_coarse    
       read(10) dummy
       read(10) omega_m,omega_l,omega_k,omega_b,dummy
       ! expansion factor, da/dtau
       read(10) aexp_ram,hexp
       close(10)
    end if

    nomfich = trim(repository)//'info_'//trim(nchar)//'.txt'
    inquire(file=nomfich,exist=ok)
    if (.not. ok) then
       write(errunit,*)'File '//trim(nomfich)//' not found'
       stop
    else
       open(unit=10,file=nomfich,status='old',form='formatted')
       do
          read (10,'(a)',end=2) line
          i = scan(line,'=')
          if (i == 0 .or. line(1:1) == '#') cycle
          name  = trim(adjustl(line(:i-1)))
          value = trim(adjustl(line(i+1:)))
          ! check for a comment at end of line !
          i     = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('unit_l' , 'scale_l')
             read(value,*) scale_l
          case ('unit_d' , 'scale_d')
             read(value,*) scale_d
          case ('unit_t' , 'scale_t')
             read(value,*) scale_t
          end select
       end do
2      close(10)
    endif

    Lboxp          = boxlen*scale_l/3.08e24/aexp_ram ! converts cgs to Mpc comoving
    aexp           = aexp_ram*af  
    omega_f        = omega_m
    omega_lambda_f = omega_l
    omega_c_f      = omega_k

    ! now read the particle data files
    nomfich = trim(repository)//'/part_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok) ! verify input file
    if ( .not. ok ) then
       write(errunit,*) trim(nomfich)//' not found.'
       stop
    endif
    open(unit=1,file=nomfich,status='old',form='unformatted')
    read(1) ncpu
    read(1) ndim
    close(1)
    
    npart = 0
    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       !write(errunit,*)'> Reading file '//trim(nomfich)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       read(1) idum
       read(1) nstar
       close(1)
       npart = npart+npart2
    end do
    
    npart   = nstar
    nbodies = npart
    write(errunit,*)'> Found ',npart,' stellar particles'
    write(errunit,*)'> Reading positions, velocities and masses...'
    
    allocate(pos(1:npart,1:ndim))
    allocate(vel(1:npart,1:ndim))
    allocate(mass(1:npart))
    ! jeje 
    allocate(ParticleID(npart))
    particleID = -1
    ! end jeje 
    
#ifdef STARS 
    allocate(stellarAge(npart),stellarZ(npart))
    stellarZ = -1.0 ! initialise metallicity to silly value 
#endif
    cnt = 0

    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)

       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       read(1) idum
       read(1) nstar
       read(1) dummy   
       read(1) dummy
       read(1) nsink

       allocate(tmpp(1:npart2,1:ndim),tmpv(1:npart2,1:ndim),tmpm(1:npart2),tmpt(1:npart2),idp(1:npart2),levelp(1:npart2))
#ifdef METALS
       allocate(dumout(1:npart2))
#endif          
       ! read all particle positions
       do idim = 1,ndim
          read(1) tmpp(1:npart2,idim)
       end do

       ! read all particle velocities
       do idim = 1,ndim
          read(1) tmpv(1:npart2,idim)
       end do

       ! read all particle masses
       read(1) tmpm(1:npart2)

       ! read all particle ids
       read(1) idp(1:npart2)

       read(1) levelp(1:npart2)

       ! read all particle creation times if necessary
       tmpt = 0.0
       if (nstar .gt. 0) then 
          read(1) tmpt(1:npart2)
#ifdef METALS
          read(1) dumout(1:npart2)
#endif          
       endif
       
       close(1)

       ! now sort particles in ascending id order and get rid of dm particles
       do ipar=1,npart2 
          if (idp(ipar) > 0 .and. tmpt(ipar) .ne. 0.0) then 
             cnt = cnt + 1
             ! put all positions between -0.5 and 0.5
             pos(idp(ipar),1:ndim) = real(tmpp(ipar,1:ndim),4) - 0.5
             ! convert code units to km/s 
             vel(idp(ipar),1:ndim) = real(tmpv(ipar,1:ndim),4)*scale_l/scale_t*1e-5
             mass(idp(ipar))       = real(tmpm(ipar),4)
             particleID(idp(ipar)) = idp(ipar)
#ifdef STARS 
             stellarAge(idp(ipar)) = tmpt(ipar)
#ifdef METALS
             stellarZ(idp(ipar))   = dumout(ipar)
#endif
#endif
          end if
       end do
       deallocate(tmpp,tmpv,tmpm,tmpt,idp,levelp)
#ifdef METALS
       deallocate(dumout)
#endif          
    end do

    if (cnt /= nstar) then 
       write(errunit,*) 'ERROR in nb of stars : (cnt/nstar)',cnt,nstar
       stop
    end if

    ! re-scale masses so that they are defined in box units (i.e. real mass = mass * mboxp)
    ! NB: for stars, mboxp is the baryonic mass of the box (Lbox^3 * rho_c * Omega_b)
    mass = mass * (scale_d * scale_l**3 *gramm_to_1011Msun) ! masses in units of 10^11Msun
    print*,'min and max masses (10^11Msun):',minval(mass),maxval(mass)
    mass  = mass / mboxp  ! where mboxp is baryonic mass of the box in 10^11Msun -> mass of particles in box mass units
    massp = maxval(mass) / nbodies  !NB: this should not be used since we leave mass allocated.

#ifdef RENORM
    write(errunit,*) 'RENORM cant be used for stars... '
    stop
#endif

!!$    open(unit=133,file='allStarDump.dat',status='unknown',form='formatted')
!!$    do ipar = 1,nbodies
!!$       write(133,'(3(e14.6,1x))') pos(ipar,1),pos(ipar,2),pos(ipar,3)
!!$    end do
!!$    close(133)

    return
    
  end subroutine read_ramses_new_stars
  
!****************************************************************************************************
  subroutine read_ramses_new(repository)
    
    ! This routine reads DM particles dumped in the RAMSES 3.0 new i/o format.  

    ! repository is directory containing output files
    ! e.g. /horizon1/teyssier/ramses_simu/boxlen100_n256/output_00001/

    implicit none

    character(len=*)            :: repository
    integer(kind=4)             :: ndim,npart,idim,icpu,ipos,ncpu,i,ipar
    integer(kind=4)             :: ncpu2,npart2,ndim2,idum,nout,nsink,nstar
    integer(kind=4)             :: nx,ny,nz,nlevelmax,ngridmax,nstep_coarse
    integer(kind=4),allocatable :: levelp(:)
    integer(kind=4),allocatable :: idp(:)
    real(kind=8)                :: boxlen,tco,aexp_ram,hexp
    real(kind=8)                :: omega_m,omega_l,omega_k,omega_b
    real(kind=8)                :: scale_l,scale_d,scale_t,dummy
    real(kind=8)                :: mtot,massres
    real(kind=8),allocatable    :: dumout(:),tmpp(:,:),tmpv(:,:),tmpm(:),tmpt(:)
    character(len=200)          :: line,name,value,nomfich
    character(len=5)            :: nchar,ncharcpu
    logical                     :: ok
    integer(kind=4) :: cnt,cntneg
    real(kind=8),parameter :: gramm_to_1011Msun = 1.d0 / 1.9891d+33 / 1.d11


    ! read cosmological params in header of amr file
    ipos    = index(repository,'output_')
    nchar   = repository(ipos+7:ipos+13)
    nomfich = trim(repository)//'amr_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok)
    if (.not. ok) then
       write(errunit,*)'File '//trim(nomfich)//' not found'
       stop
    else
       open(unit=10,file=nomfich,status='old',form='unformatted')
       read(10) ncpu
       read(10) ndim
       read(10) nx,ny,nz
       read(10) nlevelmax
       read(10) ngridmax
       read(10) idum
       read(10) idum
       read(10) boxlen
       read(10) nout,idum,idum
       allocate(dumout(nout))
       read(10) dumout
       read(10) dumout
       deallocate(dumout)
       read(10) tco
       allocate(dumout(nlevelmax))
       read(10) dumout
       read(10) dumout
       deallocate(dumout)
       read(10) idum,nstep_coarse    
       read(10) dummy
       read(10) omega_m,omega_l,omega_k,omega_b,dummy
       ! expansion factor, da/dtau
       read(10) aexp_ram,hexp
       close(10)
    end if

    nomfich = trim(repository)//'info_'//trim(nchar)//'.txt'
    inquire(file=nomfich,exist=ok)
    if (.not. ok) then
       write(errunit,*)'File '//trim(nomfich)//' not found'
       stop
    else
       open(unit=10,file=nomfich,status='old',form='formatted')
       do
          read (10,'(a)',end=2) line
          i = scan(line,'=')
          if (i == 0 .or. line(1:1) == '#') cycle
          name  = trim(adjustl(line(:i-1)))
          value = trim(adjustl(line(i+1:)))
          ! check for a comment at end of line !
          i     = scan(value,'!')
          if (i /= 0) value = trim(adjustl(value(:i-1)))
          select case (trim(name))
          case ('unit_l' , 'scale_l')
             read(value,*) scale_l
          case ('unit_d' , 'scale_d')
             read(value,*) scale_d
          case ('unit_t' , 'scale_t')
             read(value,*) scale_t
          end select
       end do
2      close(10)
    endif

    Lboxp          = boxlen*scale_l/3.08e24/aexp_ram ! converts cgs to Mpc comoving
    aexp           = aexp_ram*af  
    omega_f        = omega_m
    omega_lambda_f = omega_l
    omega_c_f      = omega_k

    ! now read the particle data files
    nomfich = trim(repository)//'/part_'//trim(nchar)//'.out00001'
    inquire(file=nomfich,exist=ok) ! verify input file
    if ( .not. ok ) then
       write(errunit,*) trim(nomfich)//' not found.'
       stop
    endif
    open(unit=1,file=nomfich,status='old',form='unformatted')
    read(1) ncpu
    read(1) ndim
    close(1)
    
    npart = 0
    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)
       !write(errunit,*)'> Reading file '//trim(nomfich)
       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       read(1) idum
       read(1) nstar
       close(1)
       npart = npart+npart2
    end do
    npart   = npart - nstar
    nbodies = npart

    write(errunit,*)'> Found ',npart,' particles'
    write(errunit,*)'> Reading positions, velocities and masses...'
    
    allocate(pos(1:npart,1:ndim))
    allocate(vel(1:npart,1:ndim))
    allocate(mass(1:npart))
    ! jeje 
    allocate(ParticleID(npart))
    particleID = -1
    ! end jeje 

    cnt = 0
    cntneg = 0
    do icpu = 1,ncpu
       call title(icpu,ncharcpu)
       nomfich = trim(repository)//'/part_'//trim(nchar)//'.out'//trim(ncharcpu)

       open(unit=1,file=nomfich,status='old',form='unformatted')
       read(1) ncpu2
       read(1) ndim2
       read(1) npart2
       read(1) idum
       read(1) nstar
       read(1) dummy   
       read(1) dummy
       read(1) nsink

       allocate(tmpp(1:npart2,1:ndim),tmpv(1:npart2,1:ndim),tmpm(1:npart2),tmpt(1:npart2),idp(1:npart2),levelp(1:npart2))
#ifdef METALS
       allocate(dumout(1:npart2))
#endif          
       ! read all particle positions
       do idim = 1,ndim
          read(1) tmpp(1:npart2,idim)
       end do

       ! read all particle velocities
       do idim = 1,ndim
          read(1) tmpv(1:npart2,idim)
       end do

       ! read all particle masses
       read(1) tmpm(1:npart2)

       ! read all particle ids
       read(1) idp(1:npart2)

       read(1) levelp(1:npart2)

       ! read all particle creation times if necessary
       tmpt = 0.0
       if (nstar .gt. 0) then 
          read(1) tmpt(1:npart2)
#ifdef METALS
          read(1) dumout(1:npart2)
#endif          
       endif
       
       close(1)

       ! now sort DM particles in ascending id order and get rid of stars !!!
       do ipar=1,npart2 
          if (idp(ipar) < 0) cntneg = cntneg + 1
          if (tmpt(ipar) == 0.0) then 
             cnt = cnt + 1
             ! put all positions between -0.5 and 0.5
             pos(idp(ipar),1:ndim) = real(tmpp(ipar,1:ndim),4) - 0.5
             ! convert code units to km/s 
             vel(idp(ipar),1:ndim) = real(tmpv(ipar,1:ndim),4)*scale_l/scale_t*1e-5
             mass(idp(ipar))       = real(tmpm(ipar),4)
             particleID(idp(ipar)) = idp(ipar)
          end if
       end do
       deallocate(tmpp,tmpv,tmpm,tmpt,idp,levelp)
#ifdef METALS
       deallocate(dumout)
#endif          
    end do
    
    ! re-allocate arrays with new nbodies = cnt
    ! ('cause perhaps there might be debris particles around, not counted as stars ...) 
    if (nbodies /= cnt) then 
       nbodies = cnt
       npart   = nbodies
       allocate(tmpp(1:nbodies,1:ndim),tmpv(1:nbodies,1:ndim),tmpm(1:nbodies),idp(1:nbodies))
       tmpp = pos(1:cnt,:)
       tmpv = vel(1:cnt,:)
       tmpm = mass(1:cnt)
       idp  = particleID(1:cnt)
       deallocate(pos,vel,mass,particleID)
       allocate(pos(1:nbodies,1:ndim),vel(1:nbodies,1:ndim),mass(1:nbodies),ParticleID(1:nbodies))
       pos        = tmpp
       vel        = tmpv
       mass       = tmpm
       particleID = idp
       deallocate(tmpp,tmpv,tmpm,idp)
    end if
    print*,nbodies,cntneg

    ! convert masses to box units (box mass is defined in compute_halo_props:init_cosmo
    mass  = mass * (scale_d * scale_l**3 * gramm_to_1011Msun)  ! DM particle mass in 10^11 Msun
    mass  = mass / mboxp ! where mboxp is total mass of the box in 10^11Msun
    massp = minval(mass)
    write(errunit,*) '> min particle mass (in 10d11 Msun)               = ',massp*mboxp
    write(errunit,*) '> max particle mass (in 10d11 Msun)               = ',maxval(mass) * mboxp

!!$    mtot = 0.0d0
!!$    do i = 1,npart
!!$       mtot = mtot + real(mass(i),8)
!!$    enddo
!!$    ! that is for the dark matter so let's add baryons now if there are any 
!!$    ! and renormalization flag is on !!
!!$    massres = minval(mass)*mboxp*1d11
!!$    massp   = minval(mass)
!!$    write(errunit,*) '> particle mass (in M_sun)               = ',massres
!!$#ifdef RENORM
!!$    massres = minval(mass)*mboxp*1d11/mtot
!!$    massp   = minval(mass)/real(mtot,4)
!!$    write(errunit,*) '> particle mass (in M_sun) after renorm  = ',massres
!!$#endif
#ifdef RENORM
    write(errunit,*) '> RENORM to be re-implemented later ... '
    stop
#endif

#ifdef BIG_RUN
    if (minval(mass) /= maxval(mass)) then 
       write(errunit,*) 'we have particles with different masses ... '
       write(errunit,*) minval(mass), maxval(mass)
       stop
    end if
    deallocate(mass)
#endif
    
  end subroutine read_ramses_new

!***********************************************************************
  subroutine write_tree_brick(nsub,nhalo,in_the_box)

    ! This subroutine writes the information relevant to building a halo 
    ! merging tree (using the build_tree program) i.e. for each halo:
    !   1/ the list of all the particles it contains (this enables us --- as  
    !      particle numbers are time independent --- to follow the halo history) 
    !   2/ its properties which are independent of merging history (mass ...)
    ! jeje : modify h%my_number just before call to write_halo, and also filter out 
    ! dubious haloes if method is not fof ... 
    ! jeje : added optional arguments ... 

    implicit none

    ! jeje 
    integer(kind=4),optional   :: nsub,nhalo
    logical(kind=4),optional   :: in_the_box(nb_of_halos+nb_of_subhalos)
    ! end jeje 
    integer(kind=4)                                         :: i,unitfile,start,j  
#ifndef BIG_RUN
    character(len=len_trim(data_dir)+16)                    :: file
#endif
    character(len=len_trim(data_dir)+len_trim(file_num)+11) :: filename
    integer(kind=4) ,allocatable                            :: members(:)
#ifdef STARS
    real(kind=4),allocatable :: tmp(:)
#endif
    logical                                                 :: done
    ! jeje 
    integer(kind=4)             :: ih,isub
!!$    integer(kind=4),allocatable :: old2new(:)
    ! end jeje 

    done = .false.
    ! leo: file now generated once with snap2mass
!!#ifndef BIG_RUN
!!    if (write_resim_masses) then 
!!       write(file,'(a,a)') trim(data_dir),'resim_masses.dat'
!!       unitfile = 44
!!       open(unitfile,file=file,form='unformatted',status='unknown')
!!       write(unitfile) nbodies
!!       write(unitfile) mass
!!       close(unitfile)     
!!       write_resim_masses = .false.
!!    end if
!!#endif

    if(.not.fsub) then
       write(filename,'(a,a,a3)') trim(data_dir),'tree_brick_',file_num
    else
       write(filename,'(a,a,a3)') trim(data_dir),'tree_bricks',file_num
    end if
    unitfile = 44
    write(errunit,*)
    write(errunit,*) '> Output data to build halo merger tree to: ',filename
    open(unit=unitfile,file=filename,form='unformatted',status='unknown')
    ! jeje : _full_sim
    write(unitfile) nbodies_full_sim
    write(unitfile) massp
    write(unitfile) aexp
    write(unitfile) omega_t
    write(unitfile) age_univ

    ! jeje : also output simply properties to a file
#ifdef SIMPLE_OUTPUTS
#ifndef STARS
    write(filename,'(a,a,a3)') trim(data_dir),'haloProps.',file_num
    open(unit=haloPropsUnit,file=filename,form='formatted',status='unknown')
#else
    write(filename,'(a,a,a3)') trim(data_dir),'galProps.',file_num
    open(unit=haloPropsUnit,file=filename,form='formatted',status='unknown')
#endif
#endif
    ! end jeje

    ! jeje : perform renumbering for sub-box (just to cope with haloes out of selected sub-box). 
    if (subboxFoF) then 
       write(unitfile) nhalo, nsub
#ifdef SIMPLE_OUTPUTS
       write(haloPropsUnit,*) nhalo,nsub
#endif
!!$       ! create index for renumbering later
!!$       allocate(old2new(nb_of_halos+nb_of_subhalos))
!!$       old2new = -1
!!$       ih   = 0
!!$       isub = nhalo
!!$       do i = 1,nb_of_halos
!!$          if (in_the_box(i)) then 
!!$             ih = ih + 1 
!!$             old2new(i) = ih 
!!$             j = liste_halos(i)%nextsub
!!$             do while (j > 0) 
!!$                isub = isub + 1
!!$                old2new(j) = isub
!!$                j = liste_halos(j)%nextsub
!!$             end do
!!$          end if
!!$       end do
    else
       write(unitfile) nb_of_halos, nb_of_subhalos
#ifdef SIMPLE_OUTPUTS
       write(haloPropsUnit,*) nb_of_halos,nb_of_subhalos
#endif
    end if
    ! end jeje

    ! jeje : output only complete halos, with continuous numbering (given by old2new array)
    if (subboxFoF) then 
       do i=1,nb_of_halos+nb_of_subhalos
          if (in_the_box(i)) then 
             allocate(members(1:nb_of_parts(i)))
             start = first_part(i)
             do j =1,nb_of_parts(i)            
                members(j) = particleID(start)
                start = linked_list(start)
             end do
             write(unitfile) nb_of_parts(i)
             write(unitfile) members
             deallocate(members)
             ! write halo properties
             liste_halos(i)%my_number = old2new(i)
             liste_halos(i)%hosthalo  = old2new(liste_halos(i)%hosthalo)
             if (liste_halos(i)%hostsub > 0) liste_halos(i)%hostsub = old2new(liste_halos(i)%hostsub)
             if (liste_halos(i)%nextsub > 0) liste_halos(i)%nextsub = old2new(liste_halos(i)%nextsub)
             call write_halo(liste_halos(i),unitfile)
#ifdef SIMPLE_OUTPUTS
             call write_halo_simple(liste_halos(i),haloPropsUnit,nb_of_parts(i))
#endif
          end if
       end do
!!$       deallocate(old2new)
    else  ! full box case 
       do i=1,nb_of_halos + nb_of_subhalos
          ! write list of particles in each halo
          allocate(members(1:nb_of_parts(i)))
          start = first_part(i)
          do j =1,nb_of_parts(i)            
             members(j) = start
             start = linked_list(start)
          end do
          write(unitfile) nb_of_parts(i)
          write(unitfile) members
          deallocate(members)
#ifdef STARS 
          ! also output masses and ages and metallicities
          allocate(tmp(1:nb_of_parts(i)))
          start = first_part(i)
          do j = 1,nb_of_parts(i)
             tmp(j) = mass(start)
             start = linked_list(start)
          end do
          write(unitfile) (tmp(j),j=1,nb_of_parts(i))
          start = first_part(i)
          do j = 1,nb_of_parts(i)
             tmp(j) = stellarAge(start)
             start = linked_list(start)
          end do
          write(unitfile) (tmp(j),j=1,nb_of_parts(i))
          start = first_part(i)
          do j = 1,nb_of_parts(i)
             tmp(j) = stellarZ(start)
             start = linked_list(start)
          end do
          write(unitfile) (tmp(j),j=1,nb_of_parts(i))
          deallocate(tmp)
#endif

          ! write each halo properties
          call write_halo(liste_halos(i),unitfile)
#ifdef SIMPLE_OUTPUTS
          call write_halo_simple(liste_halos(i),haloPropsUnit,nb_of_parts(i))
#endif
       enddo
    end if
    ! end jeje 
    close(unitfile)
#ifdef SIMPLE_OUTPUTS
    close(haloPropsUnit)
#endif

    return

  end subroutine write_tree_brick

!***********************************************************************
  subroutine write_halo(h,unitfile)

    implicit none

    integer(kind=4) :: unitfile
    type (halo)     :: h

    ! Masses (h%m,h%datas%mvir) are in units of 10^11 Msol, and 
    ! Lengths (h%p%x,h%p%y,h%p%z,h%r,h%datas%rvir) are in units of Mpc
    ! Velocities (h%v%x,h%v%y,h%v%z,h%datas%cvel) are in km/s
    ! Energies (h%ek,h%ep,h%et) are in
    ! Temperatures (h%datas%tvir) are in K
    ! Angular Momentum (h%L%x,h%L%y,h%L%z) are in
    ! Other quantities are dimensionless (h%my_number,h%my_timestep,h%spin)  

! jeje 
    write(unitfile) h%my_number!,h%minPartID
! end jeje 
    write(unitfile) h%my_timestep 
    write(unitfile) h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    write(unitfile) h%m
    write(unitfile) h%p%x,h%p%y,h%p%z
    write(unitfile) h%v%x,h%v%y,h%v%z
    write(unitfile) h%L%x,h%L%y,h%L%z 
    write(unitfile) h%r, h%sh%a, h%sh%b, h%sh%c
    write(unitfile) h%ek,h%ep,h%et
    write(unitfile) h%spin
    write(unitfile) h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    write(unitfile) h%halo_profile%rho_0,h%halo_profile%r_c

    return

  end subroutine write_halo

  !***********************************************************************
  ! jeje 
#ifdef SIMPLE_OUTPUTS
  subroutine write_halo_simple(h,unitfile,nbp)

    implicit none

    integer(kind=4) :: unitfile,nbp
    type (halo)     :: h

    write(unitfile,'(8(i16,1x),7(e14.6,1x))') h%my_number,h%minPartID,nbp,h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub,&
         & h%m,h%p%x,h%p%y,h%p%z,h%v%x,h%v%y,h%v%z

    return

  end subroutine write_halo_simple
  ! end jeje 
#endif
  !***********************************************************************
  subroutine read_last_brick(filelast) 

    implicit none
    integer(kind=4) :: idummy,nh_old,nsub_old,i,npstruct,j
    real(kind=4)    :: rdummy
    integer(kind=4),allocatable :: members(:) 

    character(len=14) ::filelast
    character(len=(len_trim(data_dir)+len_trim(filelast)))   :: filename

    ex_liste_parts = 0

    write(filename,'(a,a)') trim(data_dir),trim(filelast)
    write(errunit,*) 'Reading tree_brick from previous running:',trim(filename)
    
    open(unit=30,file=filename,form='unformatted',status='unknown')
    read(30) idummy                                         ! nbodies
    if(verbose) write(errunit,*) 'nbodies  :',idummy
    read(30) rdummy                                         ! massp
    if(verbose) write(errunit,*) 'massp    :',rdummy
    read(30) rdummy                                         ! aexp
    if(verbose) write(errunit,*) 'aexp     :',rdummy
    read(30) rdummy                                         ! omega_t
    if(verbose) write(errunit,*) 'omega_t  :',rdummy
    read(30) rdummy                                         ! age_univ
    if(verbose) write(errunit,*) 'age_univ :',rdummy
    read(30) nh_old,nsub_old                                ! nb_of_halos, nb_of_subhalos
    ex_nb_of_structs = nh_old + nsub_old
    allocate(ex_level(ex_nb_of_structs)) !! 
    do i=1,ex_nb_of_structs
       read(30)  npstruct ! nb_of_parts(i)
       !ex_nb_of_parts(i) = npstruct
       ! read list of particles in each halo
       allocate(members(npstruct))
       read(30) members
       do j =1,npstruct           
          ex_liste_parts(members(j)) = i
       end do
       deallocate(members)
       ! write each halo properties
       call read_halo(30,i)         
    enddo

    close(30)
    write(errunit,*)

    return

  end subroutine read_last_brick

  !***********************************************************************
  subroutine read_halo(unitfile,ih) 

    implicit none
    integer(kind=4) :: idummy,unitfile,ih
    real(kind=4)    :: rdummy

    read(unitfile) idummy                               !h%my_number
    read(unitfile) idummy                               !h%my_timestep 
    read(unitfile) ex_level(ih),idummy,idummy,idummy,idummy   !h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    read(unitfile) rdummy                               !h%m
    read(unitfile) rdummy,rdummy,rdummy                 !h%p%x,h%p%y,h%p%z
    read(unitfile) rdummy,rdummy,rdummy                 !h%v%x,h%v%y,h%v%z
    read(unitfile) rdummy,rdummy,rdummy                 !h%L%x,h%L%y,h%L%z 
    read(unitfile) rdummy,rdummy,rdummy,rdummy          !h%r, h%sh%a, h%sh%b, h%sh%c
    read(unitfile) rdummy,rdummy,rdummy                 !h%ek,h%ep,h%et
    read(unitfile) rdummy                               !h%spin
    read(unitfile) rdummy,rdummy,rdummy,rdummy          !h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    read(unitfile) rdummy,rdummy                        !h%halo_profile%rho_0,h%halo_profile%r_c

    return

  end subroutine read_halo

!***********************************************************************
#ifdef RECOMPUTE_PROPS
  subroutine read_tree_brick

    implicit none

    integer(kind=4)                                         :: i,unitfile,np,ip
    character(len=len_trim(data_dir)+len_trim(file_num)+11) :: filename
    integer(kind=4),allocatable                             :: plist(:)
    integer(kind=4)                                         :: prevp
    real(kind=4)                                            :: dum

    if (fsub) then
       write(filename,'(a,a,a3)') trim(data_dir),'tree_bricks',file_num
    else
       write(filename,'(a,a,a3)') trim(data_dir),'tree_brick_',file_num
    end if
    unitfile = 44
    open(unit=unitfile,file=filename,form='unformatted',status='old')
    read(unitfile) nbodies
    read(unitfile) dum ! massp
    read(unitfile) dum ! aexp
    read(unitfile) dum ! omega_t
    read(unitfile) dum ! age_univ
    read(unitfile) nb_of_halos, nb_of_subhalos

    ! initialise particle and halo lists
    allocate(linked_list(0:nbodies+1))
    allocate(first_part(0:(nb_of_halos+nb_of_subhalos)))
    allocate(nb_of_parts(0:(nb_of_halos+nb_of_subhalos)))
    allocate(liste_halos(nb_of_halos+nb_of_subhalos))
    do i = 1,nb_of_halos + nb_of_subhalos
       call clear_halo(liste_halos(i))
    end do
    first_part(0:(nb_of_halos+nb_of_subhalos))  = -1
    nb_of_parts(0:(nb_of_halos+nb_of_subhalos)) =  0
    linked_list(0:nbodies+1)                    = -1

    ! read data
    do i=1,nb_of_halos + nb_of_subhalos
       ! read list of particles in each halo                                                                    
       read(unitfile) np
       nb_of_parts(i) = np
       allocate(plist(np))
       read(unitfile) plist
       ! rebuild linked list
       first_part(i) = plist(1)
       do ip = 1,np-1
          linked_list(plist(ip)) = plist(ip+1)
       end do
       deallocate(plist)
       ! read halo props
       call really_read_halo(unitfile,liste_halos(i))
    enddo

    close(unitfile)

    return

  end subroutine read_tree_brick
#endif
!***********************************************************************

  subroutine really_read_halo(unitfile,h)

    implicit none

    type(halo)      :: h
    integer(kind=4) :: idummy,unitfile
    real(kind=4)    :: rdummy

! jeje : add minPartID
    read(unitfile) h%my_number!,h%minPartID
! end jeje
    read(unitfile) h%my_timestep
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

  end subroutine really_read_halo

!***********************************************************************
!///////////////////////////////////////////////////////////////////////

end module input_output

