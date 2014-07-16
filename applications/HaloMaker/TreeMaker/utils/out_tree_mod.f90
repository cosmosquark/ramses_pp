module out_trees

  use tree_defs

  type tree_output
     integer(kind=4) :: haloid
     integer(kind=4) :: stid
     real(kind=4)    :: m
     integer(kind=4) :: branchid
     integer(kind=4) :: sonid
     integer(kind=4) :: hostid
     integer(kind=4) :: level
     integer(kind=4) :: frag
     integer(kind=4) :: outid
  end type tree_output

  type(tree_output),allocatable :: tree_out(:)

  integer(kind=4)    :: nb_tree_plot,nb_branch,nb_struct
  integer(kind=4), allocatable :: tree_plot(:) 
  integer(kind=4),parameter    :: nb_branch_lim = 500

contains

  !**************************************************************
  subroutine search_for_halo_numbers

    use init_tree
    implicit none
    integer(kind=4) :: ilist,st,hnum,i
    
    ! search for the tree we want to output in tree_file
    n_tree_file_out = 0
    do ilist = 1, n_tree_out
       if(tree_list(ilist)%flag_out.eq.0) then
          st    = tree_list(ilist)%my_step_number
          hnum  = tree_list(ilist)%my_halo_number
          if(st.le.nsteps.and.st.ge.1) then
             call locate_int(tsno(st)%halo_number,tsno(st)%nb_of_structs,hnum,i)
          else
             i = -1
          end if
          if(i.gt.0) then
             tree_list(ilist)%flag_out   = 1
             tree_list(ilist)%my_halo_id = i
             n_tree_file_out             = n_tree_file_out + 1
             write(errunit,*) '> located halo:',&
                  tree_list(ilist)%my_step_number,tree_list(ilist)%my_halo_number
          end if
       end if
    end do

    return

  end subroutine search_for_halo_numbers

  !**************************************************************
  subroutine end_search
    
    implicit none
    integer(kind=4)    :: ilist
    character(len=200) :: wherefile
    
    if(n_tree_left.le.0) then
       write(errunit,*) '> all merger trees have been successfully outputed'
    else
       write(errunit,*) '> could not locate the following halos'
       write(errunit,*) '> among the tree_files in "tree_files.list".'
       write(errunit,*) '> h%my_timestep,h%my_number'
       do ilist = 1, n_tree_out
          if(tree_list(ilist)%flag_out.eq.0) &
               write(errunit,'(a,i3,2x,i10)') ' >',tree_list(ilist)%my_step_number,tree_list(ilist)%my_halo_number
       end do
    end if
    
    if(n_tree_out.gt.n_tree_left) then
       write(wherefile,'(a,a)') trim(dir_out),"where_haloes_are.dat"
       write(errunit,*) '> locating halos examples in full merger trees'
       write(errunit,'(a,a,a)') ' > with file "',trim(wherefile),'".'
       open(unit=37,status="unknown",form="formatted",file=wherefile)
       do ilist = 1, n_tree_out
          if(tree_list(ilist)%flag_out.eq.-1) write(37,'(2(2x,i8),5(2x,i3))') &
               tree_list(ilist)%my_tree_number,tree_list(ilist)%my_halo_number,tree_list(ilist)%my_step_number,&
               tree_list(ilist)%my_branch_number,tree_list(ilist)%my_branch_id,tree_list(ilist)%stend,tree_list(ilist)%ststart
       end do
       close(37)
    end if
    deallocate(tree_list)
    write(errunit,*) '>----------------------------------------------------------'

    return
    
  end subroutine end_search
  
  !**************************************************************
  subroutine output_selected_trees
  
    implicit none
    integer(kind=4) :: ilist
    
    call reallocate_halo_number
      
    call select_targets_in_full_trees
    ! find tree id 
    do ilist = 1, n_tree_out
       if(tree_list(ilist)%flag_out.eq.1) then
          call output_my_tree(tree_list(ilist))
       end if
    end do
    call deallocate_halo_number
    write(errunit,*)
    write(errunit,*) '>----------------------------------------------------------'
    
    
    return

  end subroutine output_selected_trees
  
  !**************************************************************
  subroutine reallocate_halo_number

    implicit none
    
    integer(kind=4) :: st
    
    do st = 1, nsteps
       if(tsno(st)%nb_of_structs.gt.0) then
          allocate(tsno(st)%halo_number(tsno(st)%nb_of_structs))
          tsno(st)%halo_number(1:tsno(st)%nb_of_structs) = -1
       end if
    end do

    return

  end subroutine reallocate_halo_number

  !**************************************************************
  subroutine deallocate_halo_number

    implicit none
    
    integer(kind=4) :: st
    
    do st = 1, nsteps
       if(tsno(st)%nb_of_structs.gt.0) then
          deallocate(tsno(st)%halo_number)
       end if
    end do

    return

  end subroutine deallocate_halo_number
   
  !**************************************************************
  subroutine select_targets_in_full_trees

    implicit none
    integer(kind=4) :: ilist,sthalo,idhalo,ifull_tree,ihost_tree

    do ilist = 1,n_tree_out
       if(tree_list(ilist)%flag_out.eq.1) then
          sthalo     = tree_list(ilist)%my_step_number
          idhalo     = tree_list(ilist)%my_halo_id
          ifull_tree = tsno(sthalo)%tree(idhalo)%my_tree_number
          data_tree(ifull_tree)%nb_target = data_tree(ifull_tree)%nb_target + 1
         
          
          ihost_tree = tsno(nsteps)%tree(ifull_tree)%hostsub
          do while(ihost_tree.gt.0) 
             ifull_tree = ihost_tree
             data_tree(ihost_tree)%nb_target = data_tree(ihost_tree)%nb_target + 1
             ihost_tree = tsno(nsteps)%tree(ihost_tree)%hostsub
          end do
          tree_list(ilist)%my_tree_id     = ifull_tree
          tree_list(ilist)%my_tree_number = data_tree(ifull_tree)%my_number
       end if
    end do
    
    return
    
  end subroutine select_targets_in_full_trees
 
  !**************************************************************
  subroutine output_my_tree(tl)
    
    implicit none
    integer(kind=4)     :: iend,iout
    type(tree_out_info) :: tl

    write(errunit,*) 
    write(errunit,*) '>----------------------------------------------------------'
    write(errunit,*) '> going to output tree for halo',tl%my_step_number,tl%my_halo_number
    write(errunit,*) '>----------------------------------------------------------'
    write(errunit,*)

    call select_elements(tl)
    if(tl%flag_out.eq.0) return
    call init_tree_out(tl,iend,iout)
    write(errunit,*) 
    write(errunit,'(a,i3,a1,i3,a1)') ' > zooming in step :[',tl%ststart,',',tl%stend,']'
    call select_structs(tl,iend,iout)
    call out_zoom_tree(tl)
  
    tl%flag_out = -1
    n_tree_left = n_tree_left - 1
    
    return
    
  end subroutine output_my_tree
  
  !**************************************************************
  subroutine select_elements(tl)
  
    implicit none
    type(tree_out_info)  :: tl 
    integer(kind=4)      :: ihalotree,isubtree,nb_tree_all,nb_tree_min,nb_branch_all,nb_branch_min,it_plot
    logical              :: limit,alliswell
    
    
    write(errunit,*) '> tree_id         :',tl%my_tree_id
    if(tl%my_tree_id.le.0) then
       write(errunit,*) '> this tree is not to be found'
       tl%flag_out = 0
       return
    end if
    write(errunit,*) '> tree_number     :',tl%my_tree_number
  
    ihalotree     = tl%my_tree_id
    if(data_tree(ihalotree)%nb_target.le.0) then
       write(errunit,*) '> Error in select_elements for tree:',ihalotree
       write(errunit,*) '> nb_target:',data_tree(ihalotree)%nb_target
       tl%flag_out = 0
       return
    end if

    if(data_tree(ihalotree)%nb_branchs.le.0) then
       write(errunit,*) '> Error in output_my_tree'
       write(errunit,*) '> no branch there..'
       tl%flag_out = 0
       return
    end if
    if(data_tree(ihalotree)%nb_structs.le.0) then
       write(errunit,*) '> Error in output_my_tree'
       write(errunit,*) '> no halo there..'
       tl%flag_out = 0
       return
    end if  


    nb_tree_all   = 1
    nb_tree_min   = 1
    nb_branch_all = data_tree(ihalotree)%nb_branchs
    nb_branch_min = data_tree(ihalotree)%nb_branchs
    isubtree      = tsno(nsteps)%tree(ihalotree)%nextsub
    do while(isubtree.gt.0)
       nb_tree_all   = nb_tree_all + 1
       nb_branch_all = nb_branch_all +  data_tree(isubtree)%nb_branchs
       if(data_tree(isubtree)%nb_target.gt.0) then
          nb_tree_min   = nb_tree_min + 1
          nb_branch_min = nb_branch_min + data_tree(isubtree)%nb_branchs
       endif
       isubtree      = tsno(nsteps)%tree(isubtree)%nextsub
    end do
    
    limit = .false.
    if(nb_branch_all.gt.nb_branch_lim) then
       if(nb_branch_min.lt.nb_branch_all) then
          limit        = .true.
          nb_tree_plot = nb_tree_min
       else
          nb_tree_plot = nb_tree_all
       endif
    else
       nb_tree_plot = nb_tree_all
    end if

    if(limit) then 
       write(errunit,*) '> the full merger tree is limited to elements containing a target'
    end if
    
    allocate(tree_plot(nb_tree_plot))
    it_plot            = 1
    tree_plot(it_plot) = ihalotree
  !  call check_tree_data(ihalotree)
    nb_branch          = data_tree(ihalotree)%nb_branchs
    nb_struct          = data_tree(ihalotree)%nb_structs
    isubtree           = tsno(nsteps)%tree(ihalotree)%nextsub
    do while(isubtree.gt.0)
       if(limit) then
          if(data_tree(isubtree)%nb_target.gt.0) then
             it_plot            = it_plot + 1
             if(it_plot.gt.nb_tree_plot) then
                write(errunit,*) '> Error in select_elements'
                write(errunit,*) '> it_plot gt nb_tree_plot:',it_plot,nb_tree_plot
                stop
             end if
             tree_plot(it_plot) = isubtree
     !        call check_tree_data(isubtree)
             nb_branch          = nb_branch + data_tree(isubtree)%nb_branchs
             nb_struct          = nb_struct + data_tree(isubtree)%nb_structs
          end if
       else
          it_plot            = it_plot + 1
          if(it_plot.gt.nb_tree_plot) then
             write(errunit,*) '> Error in select_elements'
             write(errunit,*) '> it_plot gt nb_tree_plot:',it_plot,nb_tree_plot
             stop
          end if
          tree_plot(it_plot) = isubtree
       !   call check_tree_data(isubtree)
          nb_branch          = nb_branch + data_tree(isubtree)%nb_branchs
          nb_struct          = nb_struct + data_tree(isubtree)%nb_structs
       end if
       isubtree      = tsno(nsteps)%tree(isubtree)%nextsub
    end do
    
    ! check that we have the right number of elements
    if(it_plot.ne.nb_tree_plot) then
       write(errunit,*) '> Error in select_elements'
       write(errunit,*) '> it_plot ne nb_tree_plot:',it_plot,nb_tree_plot
       stop
    end if
    ! check that the current target is there
    if(tl%my_tree_id .ne. tree_plot(1)) then
       alliswell = .false.
       ! check that tl%my_tree_id is outputted
       do it_plot = 1, nb_tree_plot
          if(tree_plot(it_plot).eq.tl%my_tree_id) alliswell = .true.
       end do
    else
       alliswell = .true.
    end if
    if(alliswell) then
       tl%my_tree_id     = tree_plot(1)
       tl%my_tree_number = data_tree(tree_plot(1))%my_number
    else
       write(errunit,*) '> Error in select elements'
       write(errunit,*) '> tree ', tl%my_tree_id,tl%my_tree_number
       write(errunit,*) '> is not in the list of elements'
       do it_plot = 1, nb_tree_plot
          write(errunit,*) '>', it_plot,tree_plot(it_plot)
       end do
       stop
    end if
    ! check that the nb of branchs is correct
    if(limit) then
       if(nb_branch.ne.nb_branch_min) then
          write(errunit,*) '> Error in select elements'
          write(errunit,*) '> tree ', tl%my_tree_id,tl%my_tree_number 
          write(errunit,*) '> nb_branch,nb_branch_min:',nb_branch,nb_branch_min
          stop
       end if
    else
       if(nb_branch.ne.nb_branch_all) then
          write(errunit,*) '> Error in select elements'
          write(errunit,*) '> tree ', tl%my_tree_id,tl%my_tree_number
          write(errunit,*) '> nb_branch,nb_branch_all:',nb_branch,nb_branch_all
          stop 
       end if
    end if

    write(errunit,*) '> in the full merger tree'
    write(errunit,*) '> nb of elements :',nb_tree_plot 
    write(errunit,*) '> nb of branchs  :',nb_branch
    write(errunit,*) '> nb of structs  :',nb_struct

!!$    integer(kind=4) :: it_plot,itree_start,lev_start
!!$    integer(kind=4) :: itree_sub,itree_host
!!$    integer(kind=4) :: nb_branch_sub,nb_branch_host,nb_branch_sys
!!$    integer(kind=4) :: nb_tree_sub,nb_tree_sys
!!$    logical         :: alliswell
!!$    
!!$    write(errunit,*) '> selecting elements'
!!$    
!!$    nb_tree_plot = 0
!!$    ! at least we plot tl%my_tree_id merger tree
!!$    itree_start  = tl%my_tree_id
!!$    do while(nb_tree_plot.le.0) 
!!$!       write(errunit,*) '> itree_start,nb_branchs',itree_start,data_tree(itree_start)%nb_branchs
!!$       ! count nb of trees  and branches we get if we add itree_start subhalos merger trees as well
!!$       lev_start      = tsno(nsteps)%tree(itree_start)%level
!!$       nb_branch_sub  = data_tree(itree_start)%nb_branchs
!!$       nb_tree_sub    = 1
!!$       itree_sub      = tsno(nsteps)%tree(itree_start)%nextsub
!!$       do while(itree_sub.gt.0) 
!!$          nb_tree_sub   = nb_tree_sub +1
!!$          nb_branch_sub = nb_branch_sub + data_tree(itree_sub)%nb_branchs
!!$!          write(errunit,*) '> itree_sub,nb_branchs',itree_sub,data_tree(itree_sub)%nb_branchs
!!$          itree_sub     = tsno(nsteps)%tree(itree_sub)%nextsub
!!$          if(itree_sub.gt.0) then
!!$             if(tsno(nsteps)%tree(itree_sub)%level.le.lev_start) itree_sub = -1
!!$          end if
!!$       end do
!!$       
!!$       ! what if we start form itree_start host, 
!!$       itree_host     =  tsno(nsteps)%tree(itree_start)%hostsub
!!$       if(itree_host.gt.0) then
!!$          ! count nb of branches we have in the host merger tree + itree_start subhalos
!!$          nb_branch_host = data_tree(itree_host)%nb_branchs + nb_branch_sub
!!$!          write(errunit,*) '> itree_host,nb_branchs',itree_host,data_tree(itree_host)%nb_branchs
!!$          ! count nb of branches we have if we also add all host subhalos merger trees
!!$          nb_branch_sys  = data_tree(itree_host)%nb_branchs
!!$          nb_tree_sys    = 1
!!$          itree_sub      = tsno(nsteps)%tree(itree_host)%nextsub
!!$          do while(itree_sub.gt.0) 
!!$             nb_tree_sys   = nb_tree_sys + 1
!!$             nb_branch_sys = nb_branch_sys + data_tree(itree_sub)%nb_branchs
!!$!             write(errunit,*) '> itree_sys,nb_branchs',itree_sub,data_tree(itree_sub)%nb_branchs
!!$             itree_sub     = tsno(nsteps)%tree(itree_sub)%nextsub
!!$             if(itree_sub.gt.0) then
!!$                if(tsno(nsteps)%tree(itree_sub)%level.le.lev_start-1)itree_sub = -1
!!$             end if
!!$          end do
!!$       else
!!$          nb_branch_host = 0
!!$          nb_branch_sys  = 0
!!$       end if
!!$!       write(errunit,*) '> itree_start,data_tree(itree_start)%nb_branchs,nb_branch_sub,nb_branch_host,nb_branch_sys'
!!$!       write(errunit,*) '>',itree_start,data_tree(itree_start)%nb_branchs,nb_branch_sub,nb_branch_host,nb_branch_sys
!!$       write(errunit,*)
!!$       if(itree_host.gt.0.and.nb_branch_host.le.nb_branch_lim) then
!!$          ! we can at least add the host merger tree
!!$          if(nb_branch_sys.gt.nb_branch_lim) then
!!$             ! we cannot add all host merger tree subhalos tree
!!$             nb_tree_plot  = nb_tree_sub + 1
!!$          else
!!$             ! we can put everything might as well start from host and check whether it is a subhalo as well and add more info
!!$             itree_start   = itree_host
!!$          end if
!!$       else
!!$          ! we cannot add the host merger tree
!!$          itree_host    = -1
!!$          if(nb_branch_sub.le.nb_branch_lim) then
!!$             ! however we can add the subhalos merger trees
!!$             nb_tree_plot = nb_tree_sub
!!$          else
!!$             ! we can only output the main merger tree
!!$             nb_tree_plot = 1
!!$          end if
!!$       end if
!!$    end do
!!$    write(errunit,*) '> nb of elements:',nb_tree_plot
!!$
!!$    allocate(tree_plot(nb_tree_plot))
!!$    tree_plot          = -1
!!$    it_plot            = 1 
!!$    nb_struct          = 0
!!$    nb_branch          = 0
!!$    if(itree_host.gt.0.and.itree_host.ne.itree_start) then
!!$       tree_plot(it_plot) = itree_host
!!$       nb_struct          = nb_struct + data_tree(itree_host)%nb_structs
!!$       nb_branch          = nb_branch + data_tree(itree_host)%nb_branchs
!!$       it_plot            = it_plot + 1
!!$       if(it_plot.gt.nb_tree_plot) then
!!$          write(errunit,*) '> Error in select_elements'
!!$          write(errunit,*) '> it_plot gt nb_tree_plot:',it_plot,nb_tree_plot
!!$          stop
!!$       end if
!!$    end if
!!$    tree_plot(it_plot) = itree_start
!!$    nb_struct          = nb_struct + data_tree(itree_start)%nb_structs
!!$    nb_branch          = nb_branch + data_tree(itree_start)%nb_branchs
!!$    if(nb_branch_sub.le.nb_branch_lim) then
!!$       lev_start    = tsno(nsteps)%tree(itree_start)%level
!!$       itree_sub    = tsno(nsteps)%tree(itree_start)%nextsub
!!$       do while(itree_sub.gt.0) 
!!$          it_plot            = it_plot +1
!!$          if(it_plot.gt.nb_tree_plot) then
!!$             write(errunit,*) '> Error in select_elements'
!!$             write(errunit,*) '> it_plot gt nb_tree_plot:',it_plot,nb_tree_plot
!!$             stop
!!$          end if
!!$          tree_plot(it_plot) = itree_sub
!!$          nb_struct          = nb_struct + data_tree(itree_sub)%nb_structs 
!!$          nb_branch          = nb_branch + data_tree(itree_sub)%nb_branchs 
!!$          
!!$          itree_sub     = tsno(nsteps)%tree(itree_sub)%nextsub
!!$          if(itree_sub.gt.0) then
!!$             if(tsno(nsteps)%tree(itree_sub)%level.le.lev_start) itree_sub = -1
!!$          end if
!!$       end do
!!$    end if
!!$    if(it_plot.ne.nb_tree_plot) then
!!$       write(errunit,*) '> Error in select_elements'
!!$       write(errunit,*) '> it_plot ne nb_tree_plot:',it_plot,nb_tree_plot
!!$       stop
!!$    end if
!!$    if(tl%my_tree_id .ne. tree_plot(1)) then
!!$       alliswell = .false.
!!$       ! check that tl%my_tree_id is outputted
!!$       do it_plot = 1, nb_tree_plot
!!$          if(tree_plot(it_plot).eq.tl%my_tree_id) alliswell = .true.
!!$       end do
!!$       if(alliswell) then
!!$          tl%my_tree_id     = tree_plot(1)
!!$          tl%my_tree_number = data_tree(tree_plot(1))%my_number
!!$       else
!!$          write(errunit,*) '> Error in select elements'
!!$          write(errunit,*) '> tree ', tl%my_tree_id,tl%my_tree_number
!!$          write(errunit,*) '> is not in the list of elements'
!!$          do it_plot = 1, nb_tree_plot
!!$            write(errunit,*) '>', it_plot,tree_plot(it_plot)
!!$          end do
!!$          stop
!!$       end if
!!$    end if
   
    return
    
  end subroutine select_elements
  
  !**************************************************************
  subroutine init_tree_out(tl,iend,iout)
    ! create 1D struct to follow up the merger tree
    
    implicit none
    integer(kind=4)     :: iend,iout
    type(tree_out_info) :: tl
    integer(kind=4)     :: itree,ibranch,istruct,istruct_old,ihalotree
    integer(kind=4)     :: st,i
    integer(kind=4)     :: mprog,nprog,ison,ihost
    logical             :: fst_time
    

    ihalotree = tl%my_tree_id
    if(tsno(nsteps)%halo_number(ihalotree) .gt. 0) then
       fst_time = .false.
       if(tsno(nsteps)%halo_number(ihalotree).ne.1) then
          write(errunit,*) '> Error in init_tree_out'
          write(errunit,*) '> ihalotree,tsno(nsteps)%halo_number(ihalotree) :',ihalotree,tsno(nsteps)%halo_number(ihalotree)
          write(errunit,*) '> ihalotree should have tsno(nsteps)%halo_number(ihalotree) eq 1'
          stop
       end if
    else
       fst_time = .true.
    end if
    
    iend    = -1
    iout    = -1
      
    allocate(tree_out(1:nb_struct))
    tree_out(1:nb_struct)%haloid   = -1
    tree_out(1:nb_struct)%stid     = -1
    tree_out(1:nb_struct)%branchid = -1
    tree_out(1:nb_struct)%m        = 0.
    tree_out(1:nb_struct)%sonid    = -1
    tree_out(1:nb_struct)%hostid   = -1
    tree_out(1:nb_struct)%level    = -1
    tree_out(1:nb_struct)%frag     = -1
    tree_out(1:nb_struct)%outid    = -1

    istruct  = 0
    ibranch  = 0
    write(errunit,*) '> tree element,nb_of_struct,nb_of_branches'
    do itree = 1, nb_tree_plot
       ibranch = ibranch + 1
       i       = tree_plot(itree)
       write(errunit,*) '> ',i,data_tree(i)%nb_structs,data_tree(i)%nb_branchs       
       st      = nsteps
       do while(i.gt.0)
          if(tsno(st)%tree(i)%my_tree_number.ne.tree_plot(itree)) then
             write(errunit,*) '> Error in init_tree_out for tree:',tl%my_tree_id,tl%my_tree_number
             write(errunit,*) '> st,i,tree_plot(itree), tsno(st)%tree(i)%my_tree_number'
             write(errunit,*) '> ',st,i,tree_plot(itree), tsno(st)%tree(i)%my_tree_number
             stop
          end if
          if(tsno(st)%tree(i)%level.le.0) then
             write(errunit,*) '> Error in init_tree_out for tree:',tl%my_tree_id,tl%my_tree_number
             write(errunit,*) '> st,i,level:',st,i, tsno(st)%tree(i)%level
             stop
          end if
          istruct     = istruct + 1
          istruct_old = tsno(st)%halo_number(i)
          if(.not.fst_time) then
             ! make sure that everything is the same as before
             if(istruct.ne.istruct_old) then
                write(errunit,*) '> Error in init_tree_out for tree:',tl%my_tree_id,tl%my_tree_number
                write(errunit,*) '> merger tree previously outputed?',.not.fst_time
                write(errunit,*) '> Error for halo:',st,i
                write(errunit,*) '> index in full tree has changed'
                write(errunit,*) '> istruct,istruct_old:',istruct,istruct_old
                stop
             end if
          else
             if(istruct_old.gt.0) then
                write(errunit,*) '> Error in init_tree_out for tree:',tl%my_tree_id,tl%my_tree_number
                write(errunit,*) '> merger tree previously outputed',.not.fst_time
                write(errunit,*) '> Error for halo:',st,i
                write(errunit,*) '> halo has already been written in a tree'
                write(errunit,*) '> istruct,istruct_old:',istruct,istruct_old
                write(errunit,*) '> st,tree_out(istruct_old)%stid:',st,tree_out(istruct_old)%stid
                write(errunit,*) '> i,tree_out(istruct_old)%haloid:',i,tree_out(istruct_old)%haloid
                stop
             end if
          end if
          if(istruct.gt.nb_struct) then
             write(errunit,*) '> Error in init_tree_out for tree:',tl%my_tree_id,tl%my_tree_number
             write(errunit,*) '> istruct,nb_struct:',istruct,nb_struct,data_tree(tree_plot(itree))%nb_structs
             write(errunit,*) '> st,i,tree_plot(itree), tsno(st)%tree(i)%my_tree_number'
             write(errunit,*) '> ',st,i,tree_plot(itree), tsno(st)%tree(i)%my_tree_number
             stop
          end if
          
          ! init all physical data
          tree_out(istruct)%stid     = st
          tree_out(istruct)%haloid   = i
          tree_out(istruct)%m        = tsno(st)%tree(i)%m
          tree_out(istruct)%branchid = ibranch 
          tree_out(istruct)%level    = tsno(st)%tree(i)%level
          tree_out(istruct)%frag     = tsno(st)%tree(i)%frag
          ! copy id in halo_number array
          tsno(st)%halo_number(i)    = istruct
          
          mprog   = tsno(st)%tree(i)%fst_prog
          if(mprog.gt.0) then
             ! go down
             i  = mprog
             st = st -1
          else
             ! go sideway
             nprog       =  tsno(st)%tree(i)%nxt_prog
             do while(nprog.le.0.and.i.gt.0)
                ! go up
                i        =  tsno(st)%tree(i)%son
                st       = st + 1
                if(i.gt.0) nprog = tsno(st)%tree(i)%nxt_prog
             end do
             if(i.gt.0) then
                i       = nprog
                ibranch = ibranch + 1
             end if
          end if
       end do
    end do
    if(istruct.ne.nb_struct) then
       write(errunit,*) '> Error in init_tree_out for tree:',tl%my_tree_id
       write(errunit,*) '> istruct,nb_struct:',istruct,nb_struct
       stop
    end if
    if(ibranch.ne.nb_branch) then
       write(errunit,*) '> Error in init_tree_out for tree:',tl%my_tree_id
       write(errunit,*) '> ibranch,nb_branch:',ibranch,nb_branch
       stop
    end if
    
    do istruct = 1, nb_struct
       st    = tree_out(istruct)%stid
       i     = tree_out(istruct)%haloid
       if(st.eq.tl%my_step_number.and.i.eq.tl%my_halo_id) then
          iout = istruct
       end if
       ison  = tsno(st)%tree(i)%son
       ihost = tsno(st)%tree(i)%hostsub
       if(ison.gt.0) then
          tree_out(istruct)%sonid = tsno(st+1)%halo_number(ison)
       else
          tree_out(istruct)%sonid = 0
       end if
       tree_out(istruct)%hostid = 0
       if(ihost.gt.0) then
          do itree = 1,nb_tree_plot
             if(tsno(st)%tree(ihost)%my_tree_number.eq.tree_plot(itree)) then
                ! host index only if it is plotted in the full merger tree
                tree_out(istruct)%hostid = tsno(st)%halo_number(ihost)
             end if
          end do
       end if
    end do

    if(iout.lt.0) then
       write(errunit,*) '> Error in init_tree_out for tree:',tl%my_tree_id
       write(errunit,*) '> could not find halo'
       write(errunit,*) '> tl%my_step_number,tl%my_halo_number,tl%my_halo_id'
       write(errunit,*) '>',tl%my_step_number,tl%my_halo_number,tl%my_halo_id
       write(errunit,*) '> should be in tree :',tsno(tl%my_step_number)%tree(tl%my_halo_id)%my_tree_number,tl%my_tree_id
       write(errunit,*) '> we were plotting trees:'
       do itree = 1, nb_tree_plot
          write(errunit,*) '> ',itree,tree_plot(itree)
       end do
       stop
    else
       tl%my_branch_number = tree_out(iout)%branchid
       iend                = iout
       if(tl%my_step_number.lt.tl%stend) then
          st = tl%my_step_number
          do while(st.ne.tl%stend)
             ison = tree_out(iend)%sonid
             if(ison.gt.0) then
                iend = ison
                st   = tree_out(iend)%stid
             else
                tl%stend = tree_out(iend)%stid
             end if
          end do
       end if
       ihost = tree_out(iend)%hostid
       if(ihost.gt.0) iend = min(iend,iout)
    end if
    deallocate(tree_plot)
    
    if(fst_time) then
       call out_full_tree(tl%my_tree_number)
       call out_full_tree_info(tl%my_tree_number)
    end if

    return

  end subroutine init_tree_out

  !**************************************************************
  subroutine check_tree_data(tree_number)
    
    
    implicit none
    integer(kind=4) :: tree_number
    integer(kind=4) :: nstructcheck,nbranchcheck
    integer(kind=4) :: st,i,mprog,nprog  
    
    write(errunit,*) '> check_tree_data :',tree_number,data_tree(tree_number)%my_number
    
    nstructcheck  = 0
    nbranchcheck  = 0
    i       = tree_number
    st      = nsteps
    do while(i.gt.0)
       if(tsno(st)%tree(i)%my_tree_number.ne.tree_number) then
          write(errunit,*) '> Error in check_tree_data for tree:',tree_number,data_tree(tree_number)%my_number
          write(errunit,*) '> st,i,tree_number, tsno(st)%tree(i)%my_tree_number'
          write(errunit,*) '> ',st,i,tree_number, tsno(st)%tree(i)%my_tree_number
          stop
       end if
       if(tsno(st)%tree(i)%level.le.0) then
          write(errunit,*) '> Error in check_tree_data for tree:',tree_number,data_tree(tree_number)%my_number
          write(errunit,*) '> st,i,level:',st,i, tsno(st)%tree(i)%level
          stop
       end if
       nstructcheck     = nstructcheck + 1
 
       mprog   = tsno(st)%tree(i)%fst_prog
       if(mprog.gt.0) then
          ! go down
          i  = mprog
          st = st -1
       else
          ! go sideway
          nprog           =  tsno(st)%tree(i)%nxt_prog
          do while(nprog.le.0.and.i.gt.0)
             ! go up
             i            =  tsno(st)%tree(i)%son
             st           = st + 1
             if(i.gt.0) nprog = tsno(st)%tree(i)%nxt_prog
          end do
          if(i.gt.0) then
             i            = nprog
             nbranchcheck = nbranchcheck + 1
          end if
       end if
       if(nstructcheck.gt.data_tree(tree_number)%nb_structs) i = -1
    end do

    if(nstructcheck.ne.data_tree(tree_number)%nb_structs.or.nbranchcheck.ne.data_tree(tree_number)%nb_branchs) then
       write(errunit,*) '>  nstructcheck,data_tree(tree_numb)%nb_structs:',nstructcheck,data_tree(tree_number)%nb_structs
       write(errunit,*) '>  nbranchcheck,data_tree(tree_numb)%nb_branchs:',nbranchcheck,data_tree(tree_number)%nb_branchs
       stop
    else
       write(errunit,*) '> all is well'
    end if

    
    return
    
  end subroutine check_tree_data

  !**************************************************************
  subroutine select_structs(tl,iend,iout)

    implicit none
    integer(kind=4)     :: iend,iout
    type(tree_out_info) :: tl
    integer(kind=4)     :: istruct,n_out
    integer(kind=4)     :: nbranch
    integer(kind=4)     :: ihoststruct,isonstruct
    integer(kind=4)     :: st,i,stson,ison,sthost,ihost
    logical             :: test_dom

    ! reinit branch id
    tree_out(1:nb_struct)%branchid = -1
     
    ! select structures between ststart and stend with halo iend and its subhalos among them
    n_out    = 0
    nbranch  = 1
    tree_out(iend)%branchid = nbranch
    do istruct = iend,nb_struct
       test_dom = tree_out(istruct)%stid.le.tl%stend.and.tree_out(istruct)%stid.ge.tl%ststart
       if(tree_out(istruct)%branchid.le.0.and.test_dom) then
          st = tree_out(istruct)%stid
          i  = tree_out(istruct)%haloid
          ! first test son
          if(tree_out(istruct)%sonid.gt.0.and.st.lt.tl%stend) then
             isonstruct = tree_out(istruct)%sonid
             if(tree_out(isonstruct)%branchid.gt.0) then
                ison       = tree_out(isonstruct)%haloid
                stson      = tree_out(isonstruct)%stid
                if(stson.ne.st+1.or.tsno(st)%tree(i)%son.ne.ison) then
                   write(errunit,*) '> Error in select_structs'
                   write(errunit,*) '> stson,st+1:',stson,st+1
                   write(errunit,*) '> ison,tsno(st)%tree(i)%son:',ison,tsno(st)%tree(i)%son
                   stop
                end if
                  
                if(tsno(stson)%tree(ison)%fst_prog.eq.i) then
                   tree_out(istruct)%branchid = tree_out(isonstruct)%branchid
                else
                   nbranch = nbranch + 1
                   tree_out(istruct)%branchid = nbranch
                end if
             end if
          else if(tree_out(istruct)%hostid.gt.0.and.st.eq.tl%stend) then
             ihoststruct = tree_out(istruct)%hostid
             if(tree_out(ihoststruct)%branchid.gt.0) then
                ihost    = tree_out(ihoststruct)%haloid 
                sthost   = tree_out(ihoststruct)%stid
                if(sthost.ne.st.or.tsno(st)%tree(i)%hostsub.ne.ihost) then
                   write(errunit,*) '> Error in select_structs'
                   write(errunit,*) '> sthost,st:',sthost,st
                   write(errunit,*) '> ihost,tsno(st)%tree(i)%hostsub:',ihost,tsno(st)%tree(i)%hostsub
                end if
                nbranch = nbranch + 1
                tree_out(istruct)%branchid = nbranch
             end if
          end if
       endif
       if(tree_out(istruct)%branchid.gt.0) then
          n_out = n_out + 1
       end if
    end do
    
    write(errunit,*) '> nb of branches  :',nbranch
    write(errunit,*) '> nb of structs   :',n_out

    call select_branches(tl,nbranch,iout)

    return

  end subroutine select_structs
  
  !**************************************************************
  subroutine select_branches(tl,nbranch,iout)
    
    implicit none
    integer(kind=4) :: nbranch,iout
    type(tree_out_info) :: tl
    integer(kind=4) :: nbranch_out,n_out,ib,ibson,istruct
    real(kind=4)    :: mbranchmin,mbranchmax,res,logres
    integer(kind=4),allocatable :: branchid(:)
    real(kind=4),allocatable    :: mbranch(:)
    
    allocate(branchid(nbranch),mbranch(nbranch))
    branchid = -1
    mbranch  = 0.

    ! compute branch mass, i.e. maximal halo mass for branchs and its secondary branches
    ! go backward to mesure secondary branch mass before main
    do istruct = nb_struct,1,-1
       ib = tree_out(istruct)%branchid
       if(ib.gt.0) then
          mbranch(ib) = max(mbranch(ib),tree_out(istruct)%m)
       end if
    end do
    mbranchmax = maxval(mbranch)
    mbranchmin = minval(mbranch)
    

    ! give mbrnchmax mass to all branches connecting halo iout to step stend
    istruct  = iout
    ib       = -1
    ibson    = tree_out(istruct)%branchid
    do while(istruct.gt.0)
       if(ibson .ne. ib) then
          mbranch(ibson) = mbranchmax
          ib             = ibson
       end if
       istruct     = tree_out(istruct)%sonid
       if(istruct.gt.0) then
          ibson = tree_out(istruct)%branchid
          if(ibson.le.0) istruct = -1
       end if   
    end do
    mbranchmin = minval(mbranch)
    
    if(res_force.gt.0.) then
       mbranchmin  = mbranchmax*res_force
       ! always keep he main branch
       nbranch_out = 1
       branchid(1) = nbranch_out
       do ib = 2,nbranch
          if(mbranch(ib).gt.mbranchmin) then
             nbranch_out = nbranch_out + 1
             branchid(ib) = nbranch_out
          else
             branchid(ib) = -1
          end if
       end do
       res    = res_force
    else
       logres      = - 5.
       res         = 0
       nbranch_out = nbranch
       do while(nbranch_out.gt.nbranchmax)
          do while(mbranchmax*res.le.mbranchmin)
             logres      = logres + 0.25 
             res         = 10.**logres
          end do
          mbranchmin  = mbranchmax*res
          ! always keep he main branch
          nbranch_out = 1
          branchid(1) = nbranch_out
          do ib = 2,nbranch
             if(mbranch(ib).gt.mbranchmin) then
                nbranch_out = nbranch_out + 1
                branchid(ib) = nbranch_out
             else
                branchid(ib) = -1
             end if
          end do
       end do
    end if
   
    write(errunit,*) '> resolution      :',res
    write(errunit,*) '> max branch mass :',mbranchmax*1e11
    write(errunit,*) '> min branch mass :',mbranchmin*1e11
    write(errunit,*) '> nb of branchs   :',nbranch_out

    n_out = 0
    do istruct = 1,nb_struct
       ib = tree_out(istruct)%branchid
       if(ib.gt.0) then
          if(nbranch.ne.nbranch_out) then
             tree_out(istruct)%branchid = branchid(ib)
          end if
          if(tree_out(istruct)%branchid.gt.0) then
             n_out                      = n_out + 1
             tree_out(istruct)%outid    = n_out
          end if
       end if
    end do

    tl%my_branch_id = tree_out(iout)%branchid
   
    write(errunit,*) '> nb of structs   :',n_out
    deallocate(branchid,mbranch)

    return
    
  end subroutine select_branches
  
  !**************************************************************
  subroutine out_full_tree(tree_number)
  
    implicit none
    integer(kind=4)           :: tree_number
    integer(kind=4)           :: istruct
    character(len=200)        :: tree_out_file 

    write(errunit,*)
    write(errunit,*) '> outputing full merger tree'
    write(tree_out_file,'(a,a,i8.8,a1,a3,a)') &
         trim(dir_out),"full_tree_",tree_number,"_",method,".dat"
    write(errunit,'(a,a,a)') ' > in file "',trim(tree_out_file),'".'
 
    open(unit=53,form='formatted',status='unknown',file=tree_out_file)
    write(53,'(4(1x,i6),1x,E10.4,2(1x,i2))') 0,0,0,0,0.0,0,0
    do istruct = 1,nb_struct
       write(53,'(4(1x,i6),1x,E10.4,2(1x,i2))') &
            tree_out(istruct)%branchid,tree_out(istruct)%stid,tree_out(istruct)%sonid,tree_out(istruct)%hostid,&
            tree_out(istruct)%m*1e11,tree_out(istruct)%level,tree_out(istruct)%frag
    end do
    close(53)
       
    return
    
  end subroutine out_full_tree

  !**************************************************************
  subroutine out_full_tree_info(tree_number)
    
    implicit none 
    integer(kind=4)           :: tree_number
    integer(kind=4)           :: istruct,stid,haloid,halo_number,sonid,son_number
    character(len=200)        :: tree_out_file 
    
    write(errunit,*)
    write(errunit,*) '> outputing full merger tree info'
    write(tree_out_file,'(a,a,i8.8,a1,a3,a)') &
         trim(dir_out),"merger_tree_info_",tree_number,"_",method,".dat"
    write(errunit,'(a,a,a)') ' > in file "',trim(tree_out_file),'".'
    
    open(unit=63,form='formatted',status='unknown',file=tree_out_file)
    write(63,*) "branch,step,halo,son,level,frag"
    do istruct = 1,nb_struct
       stid          = tree_out(istruct)%stid
       haloid        = tree_out(istruct)%haloid
       halo_number   = tsno(stid)%tree(haloid)%my_number
       son_number    = 0
       sonid         = tree_out(istruct)%sonid
       if(sonid.gt.0) then
          sonid      = tree_out(sonid)%haloid
          son_number = tsno(stid+1)%tree(sonid)%my_number
       end if
       write(63,'(1x,i4,1x,i3,2(1x,i8.8),2(1x,i2))') &
            tree_out(istruct)%branchid,tree_out(istruct)%stid,halo_number,son_number,tree_out(istruct)%level,tree_out(istruct)%frag
    end do
    close(63)
       
    return
    
  end subroutine out_full_tree_info

  !**************************************************************
  subroutine out_zoom_tree(tl)
  
    implicit none
    type(tree_out_info)       :: tl
    integer(kind=4)           :: istruct,ison,ihost  
    character(len=200)        :: tree_out_file

    write(errunit,*)
    write(errunit,*) '> outputing zoom merger tree'
    write(tree_out_file,'(a,a,2(i8.8,a1),3(i3.3,a1),a3,a)') &
         trim(dir_out),"zoom_tree_",tl%my_tree_number,"_",     &
         tl%my_halo_number,"_",tl%my_step_number,"_",tl%ststart,"_",tl%stend,&
         "_",method,".dat"
    write(errunit,'(a,a,a)') ' > in file "',trim(tree_out_file),'".'
    
    open(unit=54,form='formatted',status='unknown',file=tree_out_file)
    write(54,'(4(1x,i6),1x,E10.4,2(1x,i2))') 0,0,0,0,0.0,0,0
    do istruct = 1,nb_struct
       if(tree_out(istruct)%outid.gt.0) then
          ison  = tree_out(istruct)%sonid
          if(ison.gt.0) ison = tree_out(ison)%outid
          ihost = tree_out(istruct)%hostid
          if(ihost.gt.0) ihost = tree_out(ihost)%outid
          write(54,'(4(1x,i6),1x,E10.4,2(1x,i2))') tree_out(istruct)%branchid,tree_out(istruct)%stid,ison,ihost,&
               tree_out(istruct)%m*1e11,tree_out(istruct)%level,tree_out(istruct)%frag
       end if
    end do
    close(54)
    
    deallocate(tree_out)
    
    return

  end subroutine out_zoom_tree
  
  !**************************************************************

end module out_trees
