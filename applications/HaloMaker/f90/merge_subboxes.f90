program merge_subboxes
  
  use input_output
  
  implicit none
  
  character(200)              :: rundir,filename
  integer(kind=4)             :: nsnaps,nsubbox,i,snap,j,k,ip,np,unitfile,isub,l,m
  integer(kind=4),allocatable :: snaps(:),plist(:),old(:)

  logical(kind=4)             :: ex,ok
  integer(kind=4)             :: ns,nh,nstot,nhtot,minid
  real(kind=4)                :: rdum

  call getarg(1,rundir)  ! directory where the outputs of HaloMaker are

  ! read list of processed snapshots (and number of sub-boxes) to assemble
  write(filename,'(a,a)') trim(rundir), '/mergeFile.sh'
  open(unit=12,file=filename,status='old',form='formatted')
  read(12,*) nsnaps,nsubbox
  allocate(snaps(nsnaps)) 
  do i = 1,nsnaps
     read(12,*) snaps(i)
  end do
  close(12)

  ! loop on snapshots and assemble them one by one .. 
  do i = 1,nsnaps
     snap = snaps(i)
     ! count total number of halos and sub-halos
     nstot = 0
     nhtot = 0
     nbodies_full_sim = 0
     do j = 0,nsubbox-1
        write(filename,'(a,a,a,a,a,a,i3.3)') trim(rundir),'/',trim(int2string(snap)),'/',trim(int2string(j)),'/tree_brick_',snap
        if (j == 0) then 
           inquire(file=filename,exist=ex)
           if (ex) then
              write(errunit,*) '> No sub-structures ... '
           else
              write(filename,'(a,a,a,a,a,a,i3.3)') trim(rundir),'/',trim(int2string(snap)),'/',trim(int2string(j)),'/tree_bricks',snap
              write(errunit,*) '> Sub-structures ... '
           end if
        else 
           if (.not. ex) write(filename,'(a,a,a,a,a,a,i3.3)') trim(rundir),'/',trim(int2string(snap)),'/',trim(int2string(j)),'/tree_bricks',snap
        end if

        ! file may not exist if there are no structures ... 
        inquire(file=filename,exist=ok)        
        if (ok) then 
           open(unit=12,file=filename,status='old',form='unformatted')
           read(12) nbodies ! of full simulation
           read(12) rdum
           read(12) rdum
           read(12) rdum
           read(12) rdum
           read(12) nh,ns
           close(12)
           nbodies_full_sim = nbodies
        else
           nh      = 0
           ns      = 0
        end if
        nhtot = nhtot + nh
        nstot = nstot + ns
     end do
     print*, 'Total number of Halos and sub-halos : ',nhtot,nstot
     nb_of_halos      = nhtot
     nb_of_subhalos   = nstot
     
     ! allocate/initialize stuff
     allocate(liste_halos(nb_of_halos+nb_of_subhalos))
     do k = 1,nb_of_halos+nb_of_subhalos
        call clear_halo(liste_halos(k))
     end do
     allocate(linked_list(0:nbodies+1))
     allocate(first_part(0:(nb_of_halos+nb_of_subhalos)))
     allocate(nb_of_parts(0:(nb_of_halos+nb_of_subhalos)))
     first_part(:)  = -1
     nb_of_parts(:) =  0
     linked_list(:) = -1
     
     ! read data
     unitfile = 44
     nhtot = 0
     nstot = nb_of_halos
     do j=0,nsubbox-1
        if (ex) then 
           write(filename,'(a,a,a,a,a,a,i3.3)') trim(rundir),'/',trim(int2string(snap)),'/',trim(int2string(j)),'/tree_brick_',snap
        else
           write(filename,'(a,a,a,a,a,a,i3.3)') trim(rundir),'/',trim(int2string(snap)),'/',trim(int2string(j)),'/tree_bricks',snap
        end if
        inquire(file=filename,exist=ok)
        if (ok) then 
           open(unit=unitfile,file=filename,status='old',form='unformatted')
           read(unitfile) nbodies
           read(unitfile) massp
           read(unitfile) aexp
           read(unitfile) omega_t
           read(unitfile) age_univ
           read(unitfile) nh,ns
           print*,nh,ns
           allocate(old2new(nh+ns))
           old2new = -1
           allocate(old(nh+ns))
           ! read main halos first
           do k=1,nh
              nhtot = nhtot + 1
              old2new(k) = nhtot
              read(unitfile) np
              nb_of_parts(nhtot) = np
              allocate(plist(np))
              read(unitfile) plist
              ! rebuild linked list
              first_part(nhtot) = plist(1)
              do ip = 1,np-1
                 linked_list(plist(ip)) = plist(ip+1)
              end do
              deallocate(plist)
              ! read halo props
              call really_read_halo(unitfile,liste_halos(nhtot))
              old(k) = liste_halos(nhtot)%my_number
           end do
           ! read sub-halos
           do k=1,ns
              nstot = nstot + 1
              old2new(k+nh) = nstot
              read(unitfile) np
              nb_of_parts(nstot) = np
              allocate(plist(np))
              read(unitfile) plist
              ! rebuild linked list
              first_part(nstot) = plist(1)
              do ip = 1,np-1
                 linked_list(plist(ip)) = plist(ip+1)
              end do
              deallocate(plist)
              ! read halo props
              call really_read_halo(unitfile,liste_halos(nstot))
              old(k+nh) = liste_halos(nstot)%my_number
           end do
           
           ! global renumbering of sub-box (sub-)halos
           do k = 1,nh+ns
              l = old2new(k)
              liste_halos(l)%my_number = l
              call locate_int(old,nh+ns,liste_halos(l)%hosthalo,m)
              !print*,liste_halos(l)%hosthalo,m,old(m),old2new(m)
              if (old(m) /= liste_halos(l)%hosthalo) then 
                 print*,'bloody hell, this is now allowed ... ',liste_halos(l)%hosthalo
                 !print*,old
                 stop
              end if
              liste_halos(l)%hosthalo  = old2new(m)
              
              if (liste_halos(l)%hostsub > 0) then 
                 call locate_int(old,nh+ns,liste_halos(l)%hostsub,m)
                 liste_halos(l)%hostsub = old2new(m)
                 !!liste_halos(l)%hostsub = old2new(liste_halos(l)%hostsub)
              end if
              if (liste_halos(l)%nextsub > 0) then 
                 call locate_int(old,nh+ns,liste_halos(l)%nextsub,m)
                 liste_halos(l)%nextsub = old2new(m)
                 !!liste_halos(l)%nextsub = old2new(liste_halos(l)%nextsub)
              end if
           end do
           deallocate(old2new,old)
           close(unitfile)
        end if
     end do ! loop on subboxes

#ifdef SIMPLE_OUTPUTS
     ! Compute minPartID
     do k = 1,nb_of_halos+nb_of_subhalos
        ip = first_part(k)
        liste_halos(k)%minPartID = ip
        do while (ip > 0) 
           if (liste_halos(k)%minPartID > ip) liste_halos(k)%minPartID = ip
           ip = linked_list(ip)
        end do
     end do
#endif

     ! re-write complete re-numbered catalog (new tree_brick)
     write(data_dir,'(a,a,a,a)') trim(rundir),'/',trim(int2string(snap)),'/'
     if (ex) then 
        fsub = .false.
     else
        fsub = .true.
     end if
     subboxFoF      = .false.
     write(file_num,'(i3.3)') snap
     call write_tree_brick
   
     deallocate(linked_list,first_part,nb_of_parts)
     deallocate(liste_halos)

  end do ! loop on snapshots
  
  deallocate(snaps)

contains

  function int2string(i)
    implicit none 
    integer(kind=4) :: i
    character(4)   :: int2string
    if (i<10) then 
       write(int2string,'(i1)') i
    else if (i < 100) then 
       write(int2string,'(i2)') i
    else if (i < 1000) then 
       write(int2string,'(i3)') i
    else if (i < 10000) then 
       write(int2string,'(i4)') i
    else 
       stop 'Hey!'
    end if
    return
  end function int2string

  subroutine print_halo(h)
    implicit none
    type(halo) :: h 
    print*,h%my_number
    print*,h%my_timestep 
    print*,h%level,h%hosthalo,h%hostsub,h%nbsub,h%nextsub
    print*,h%m
    print*,h%p%x,h%p%y,h%p%z
    print*,h%v%x,h%v%y,h%v%z
    print*,h%L%x,h%L%y,h%L%z 
    print*,h%r, h%sh%a, h%sh%b, h%sh%c
    print*,h%ek,h%ep,h%et
    print*,h%spin
    print*,h%datas%rvir,h%datas%mvir,h%datas%tvir,h%datas%cvel
    print*,h%halo_profile%rho_0,h%halo_profile%r_c
    return
  end subroutine print_halo

  !*****************************************************************************************************************

  subroutine locate_int(xx,n,x,j)

    ! subroutine which locates the position of a value x in an array xx chap 3

    implicit none

    integer(kind=4) ::  n,j,jl,ju,jm
    integer(kind=4) ::  xx(n),x

    do j = 1,n
       if (xx(j) == x) exit
    end do
    
!!$    jl = 1
!!$    ju = n+1
!!$    do while (ju-jl > 1) 
!!$       jm = (ju+jl)/2
!!$       if (x > xx(jm)) then
!!$          jl = jm
!!$       else
!!$          ju = jm
!!$       endif
!!$    enddo
!!$    j = jl

    return

  end subroutine locate_int

end program merge_subboxes


