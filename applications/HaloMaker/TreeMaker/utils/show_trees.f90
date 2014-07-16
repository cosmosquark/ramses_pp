program show_my_trees

  use tree_defs
  use input_output
  use init_tree
  use out_trees
  
  implicit none
  integer(kind=4)    :: ifile,n_tree_files_in,ierr_file_list
  character(len=200) :: file_list,file_tree_loc
  
  write(errunit,*) '>===================================='
  write(errunit,*) 
  write(errunit,*) ' Show Tree '
  write(errunit,*) 
  write(errunit,*) '>===================================='
  write(errunit,*) 
   
  call getarg(1,data_dir)       ! get directory where to input/output the data
  call read_init

  write(file_list,'(a,a)') trim(dir_tree),"tree_files.list" 
  open(unit=72,form="formatted",status="old",file=file_list,iostat=ierr_file_list)
  if(ierr_file_list.ne.0) then
     write(errunit,*) '> file "',trim(file_list),'" doesn''t exist'
     ! if file is not there create it
     open(unit=73,form="formatted",status="unknown",file=file_list)
     ! if no tree_file are going to be red, read first one
     if(n_tree_files.le.0) n_tree_files = 1
     write(73,*) n_tree_files
  else
     read(72,*) n_tree_files_in
     ! if no tree_file are going to be red, read all
     if(n_tree_files.le.0.or.n_tree_files.gt.n_tree_files_in) &
          n_tree_files = n_tree_files_in
  end if
  if(n_tree_files.le.0) then
     write(errunit,*) '> no files to read...'
     stop
  endif
  if(n_tree_out.le.0) then
     call open_tree_data_info
  else
     n_tree_left = n_tree_out
  end if
  do ifile = 1, n_tree_files
     if(ierr_file_list.ne.0) then
        ! no tree_file in list make a list starting from the first tree file 
        write(file_tree_loc,'(a,i3.3,a1,i3.3)') 'tree_file_',nsteps,'.',ifile
        write(73,*) trim(file_tree_loc)
     else
        read(72,*) file_tree_loc
     end if
     write(file_tree,*) trim(dir_tree),trim(file_tree_loc)
     if(n_tree_out.le.0.or.n_tree_left.gt.0) then
        call read_tree
        if(n_tree_left.gt.0) then
           call search_for_halo_numbers
        else
           n_tree_file_out = 0
        end if
        if(n_tree_file_out.gt.0.or.n_tree_out.le.0) then
           call renumber_all
           call init_tree_datas
           if(n_tree_out.le.0) then
              call write_tree_data_info
           else
              call output_selected_trees
           end if
           call deallocate_all
        else
           call deallocate_tree_only
        end if
     end if
  end do
  if(ierr_file_list.ne.0) then
     close(73)
  else
     close(72)
  end if
  
  !! make seperate routine
  !!! remember the where it is stuff.
  if(n_tree_out.gt.0) then
     call end_search
  end if
  write(errunit,*) 
  write(errunit,*) '>===================================='
  write(errunit,*) ' End '
  write(errunit,*) '>===================================='

  stop

end program show_my_trees
