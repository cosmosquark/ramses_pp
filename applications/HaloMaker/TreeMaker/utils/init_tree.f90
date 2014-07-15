module init_tree

  use tree_defs
  
contains
  
  !**************************************************************
  subroutine init_tsno
    
    implicit none
    integer(kind=4) :: st
    
    do st = 1, nsteps
       
       tsno(st)%aexp           = 0.
       tsno(st)%omega_t        = 0.
       tsno(st)%age_univ       = 0.
       tsno(st)%redshift       = 0.
       tsno(st)%nb_of_halos    = 0
       tsno(st)%nb_of_subhalos = 0
       tsno(st)%nb_of_structs  = 0
       
       nullify(tsno(st)%halo_number)
       nullify(tsno(st)%tree)
       
    end do
        
    return

  end subroutine init_tsno
  
  !**************************************************************
  subroutine init_halotree(ht)
    
    implicit none 
    type(tree_info) :: ht
    
    ht%my_number      = -1
    ht%my_tree_number = -1
    ht%m              = 0.
    ht%mvir           = 0.
    ht%level          = 0
    ht%hosthalo       = -1
    ht%hostsub        = -1
    ht%nextsub        = -1
    ht%son            = -1
    ht%ndads          = 0
    ht%fst_prog       = -1
    ht%nxt_prog       = -1
    ht%frag           = 0
    nullify(ht%dads)
    nullify(ht%mass_dads)
    
    return
    
  end subroutine init_halotree

  !**************************************************************
  subroutine init_tree_list
    
    implicit none
    
    tree_list(1:n_tree_out)%my_halo_number = -1
    tree_list(1:n_tree_out)%my_step_number = -1
    tree_list(1:n_tree_out)%my_tree_number = -1
    tree_list(1:n_tree_out)%stend          = -1
    tree_list(1:n_tree_out)%ststart        = -1
    tree_list(1:n_tree_out)%my_tree_id     = -1
    tree_list(1:n_tree_out)%my_halo_id     = -1
    tree_list(1:n_tree_out)%my_branch_id   = -1
    tree_list(1:n_tree_out)%flag_out       = 0
    
    return

  end subroutine init_tree_list

  !**************************************************************
  subroutine renumber_all

    implicit none 
    integer(kind=4) :: st,i
    type(tree_info), pointer :: ht
    

    write(errunit,*) '> renumbering...'
    do st = 1,nsteps                                ! loop on all timesteps
       if(tsno(st)%nb_of_structs.gt.0) then 
          do i = 1,tsno(st)%nb_of_structs              ! loop on all halos in a given timestep
             ht => tsno(st)%tree(i)                    ! h is an alias for current halo
             call renumber_tree(ht,st)   
             if(st.gt.1) then 
                call make_prog_linked_list(ht,st,i)
             end if
          end do
          if(st.gt.1) then
             ! no nead for tsno(st)%halo_number anymore
             if(tsno(st-1)%nb_of_structs.gt.0) &
                  deallocate(tsno(st-1)%halo_number)
          end if
       end if
    end do
    ! no nead for tsno(nsteps)%halo_number anymore
    if(tsno(nsteps)%nb_of_structs.gt.0) &
         deallocate(tsno(nsteps)%halo_number)
    
    return

  end subroutine renumber_all

 !***************************************************************
  subroutine renumber_tree(ht,st)
    
    implicit none
    integer(kind=4) :: st
    type(tree_info) :: ht
    integer(kind=4) :: idad,renum
    
    if(ht%ndads.gt.0) then
       ! since we are renumbering add something to the frag flag
       if(ht%frag.ne.0) ht%frag = 2
       if(st.le.1) then
          write(errunit,*) '> Error in renumber_tree'
          write(errunit,*) '> cannot renum dads at step:',st-1
          stop
       end if
       do idad = 1,ht%ndads
          if(ht%dads(idad).le.0) then
             write(errunit,*) '> Error in renumber_tree'
             write(errunit,*) '> idad,t%dads(idad):',idad,ht%dads(idad)
             write(errunit,*) '> background should not be amoung the dads'
             stop
          end if
          renum = ht%dads(idad)
          call renumbering(renum,st-1)
          if(renum.le.0) then
             write(errunit,*) '> Error in renumber_tree'
             write(errunit,*) '> for dad:',st-1,idad,ht%dads(idad)
             stop
          else
             ht%dads(idad) = renum
          end if
       end do
    end if

    if(ht%son.gt.0) then
       if(st.ge.nsteps) then
          write(errunit,*) '> Error in renumber_tree'
          write(errunit,*) '> cannot renum sons at step:',st+1
          stop
       end if
       renum = ht%son
       if(st.lt.nsteps) then
          call renumbering(renum,st+1)
       else
          ! in case we do not use the whole simulation reset %son at steps nsteps_do
          renum = 0
       end if
       ht%son = renum
    end if

    if(ht%hosthalo.gt.0) then
       renum = ht%hosthalo
       call renumbering(renum,st)
       if(renum.le.0) then
          write(errunit,*) '> Error in renumber_tree'
          write(errunit,*) '> for hosthalo:',st,ht%hosthalo
          stop
       else
          ht%hosthalo = renum
       end if
    end if
    
    if(ht%hostsub.gt.0) then
       renum = ht%hostsub
       call renumbering(renum,st)
       if(renum.le.0) then
          write(errunit,*) '> Error in renumber_tree'
          write(errunit,*) '> for hostsub:',st,ht%hostsub
          stop
       else
          ht%hostsub = renum
       end if
    end if
    
    if(ht%nextsub.gt.0) then
       renum = ht%nextsub
       call renumbering(renum,st)
       if(renum.le.0) then
          write(errunit,*) '> Error in renumber_tree'
          write(errunit,*) '> for nextsub:',st,ht%nextsub
          stop
       else
          ht%nextsub = renum
       end if
    end if

    return
    
  end subroutine renumber_tree

  !**************************************************************
  subroutine renumbering(renum,tsno_ts)

    ! this routine reorders halos so that renum matches the tsno 
    ! list number (needed if a big simulation is split in chunks)

    implicit none 

    integer(kind=4) :: renum,tsno_ts
    integer(kind=4) :: index, n

    index = -1
    if(renum.gt.0)then
       n = tsno(tsno_ts)%nb_of_structs
       call locate_int(tsno(tsno_ts)%halo_number,n,renum,index)
       if(.not.(renum.eq.tsno(tsno_ts)%halo_number(index)))then
          write(errunit,*) '> Fatal error in renum routine'
          write(errunit,*) '> renum,tsno(tsno_ts)%tree(index)%my_number:',renum,tsno(tsno_ts)%tree(index)%my_number
          write(errunit,*) '> tsno_ts,tsno(tsno_ts)%nb_of_structs',tsno_ts,tsno(tsno_ts)%nb_of_structs     
          if(index .gt. 1 ) write(errunit,*) '> idex -1,tsno(tsno_ts)%halo_number(index-1)',index-1,tsno(tsno_ts)%halo_number(index-1)
          write(errunit,*) '> index,tsno(tsno_ts)%halo_number(index)', index,tsno(tsno_ts)%halo_number(index)
          if(index .lt. tsno(tsno_ts)%nb_of_structs ) write(errunit,*)'> index +1,tsno(tsno_ts)%halo_number(index+1)',index+1,tsno(tsno_ts)%halo_number(index+1)
          stop
       end if
       renum = index
    endif

    return

  end subroutine renumbering

  !*****************************************************************************************************************
  subroutine locate_int(xx,n,x,j)

    ! subroutine which locates the position of a value x in an array xx chap 3

    implicit none

    integer(kind=4) ::  n,j,jl,ju,jm
    integer(kind=4) ::  xx(n),x

    ! DT modification of locate_int, ass xx is ordered for x to be found in xx
    ! x .ge.xx(1) and x.le.xx(n)
    if(x.lt.xx(1).or.x.gt.xx(n)) then
       j = -1
       return
    end if
    
    jl = 1
    ju = n+1
    
    do while (ju-jl > 1) 
       jm = (ju+jl)/2
       if (x > xx(jm)) then
          jl = jm
       else
          ju = jm
       endif
    enddo
    !! DT modification of locate_int old version doesn't work with n = 3 sometimes
    if(xx(jl) .eq. x) then
       j = jl
    else if(xx(ju) .eq. x) then
       j = ju
    else
       !! DT put -1 if not found, sure to make the code stop
       j = -1
    end if
     
    return

  end subroutine locate_int

  !*****************************************************************************************************************
  subroutine make_prog_linked_list(ht,st,i)
    
    implicit none
    type(tree_info) :: ht
    integer(kind=4) :: st,i
    type(tree_info),pointer :: dad
    integer(kind=4) :: idad,iddad

    if(ht%ndads.le.0) return
    
    if(ht%ndads.eq.1) then
       ht%fst_prog = ht%dads(1)
    else
       do idad  = ht%ndads,1,-1
          iddad = ht%dads(idad)
          dad   => tsno(st-1)%tree(iddad)
          if(dad%son.ne.i) then
             write(errunit,*) '> Error in make_prog_list for halo',st,ht%my_number
             write(errunit,*) '> idad,iddad,dad%son,i'
             write(errunit,*) '>',idad,iddad,dad%son,i
             stop
          end if
          dad%nxt_prog = ht%fst_prog
          ht%fst_prog  = iddad 
       end do
    end if
    
    deallocate(ht%dads)
    deallocate(ht%mass_dads)
    
    return
    
  end subroutine make_prog_linked_list

  !**************************************************************
  subroutine init_tree_datas

    implicit none
    integer(kind=4) :: tree_numb,st,i,ison
    type(tree_info), pointer :: ht,son
    
    write(errunit,*) '> computing tree data...'
    
    if(tsno(nsteps)%nb_of_halos.le.0) then
       write(errunit,*) '> Error in compute_tree_datas'
       write(errunit,*) '> tsno(nsteps)%nb_of_halos:',tsno(nsteps)%nb_of_halos
       stop
    end if

    allocate(data_tree(tsno(nsteps)%nb_of_structs))
    
    write(errunit,*) '> initiate tree data at last step'
    do i = 1, tsno(nsteps)%nb_of_structs
       tree_numb = i
       ht => tsno(nsteps)%tree(i)
       data_tree(tree_numb)%my_number    = ht%my_number
       ht%my_tree_number                 = tree_numb
       data_tree(tree_numb)%nb_branchs   = 1
       data_tree(tree_numb)%nb_structs   = 1
       if(ht%frag.eq.1) then
          data_tree(tree_numb)%nb_frags  = 1
       else
          data_tree(tree_numb)%nb_frags  = 0
       end if
       data_tree(i)%nb_target            = 0
    end do

    write(errunit,*) '> computing tree_index and branch index at all step'
    do st = nsteps-1, 1, -1
       do i = 1,tsno(st)%nb_of_structs
          ht => tsno(st)%tree(i)
          ison = ht%son
          if(ison.gt.0) then
             son => tsno(st+1)%tree(ison)
             tree_numb = son%my_tree_number
             if(tree_numb.gt.0) then
                data_tree(tree_numb)%nb_structs = data_tree(tree_numb)%nb_structs + 1
                ht%my_tree_number               = tree_numb
                ! branch numbering
                if(son%fst_prog.ne.i) then
                   data_tree(tree_numb)%nb_branchs   = data_tree(tree_numb)%nb_branchs + 1
                end if
                if(ht%frag.eq.1) then
                   data_tree(tree_numb)%nb_frags  = data_tree(tree_numb)%nb_frags + 1
                end if
             end if
          end if
       end do
    end do
    
    write(errunit,*) '> checking results'
    do tree_numb = 1,tsno(nsteps)%nb_of_structs
       call check_numbers(tree_numb)
    end do
    write(errunit,*) '> So far everything is fine'
    write(errunit,*)
    
    return

  end subroutine init_tree_datas

  !*****************************************************************************************************************
  subroutine check_numbers(itree)
    
    implicit none
    integer(kind=4) :: itree
    integer(kind=4) :: st,i,nprog,mprog
    integer(kind=4) :: nstruct_ch,nbranch_ch
   
    nstruct_ch = 0
    nbranch_ch = 1
    
    st = nsteps
    i  = itree
    do while(i.gt.0)
       if(tsno(st)%tree(i)%my_tree_number.ne.itree) then
          write(errunit,*) '> Error in check_numbers for tree:',itree,data_tree(itree)%my_number
          write(errunit,*) '> st,ih,itree,tsno(st)%tree(i)%my_tree_number'
          write(errunit,*) '> ',st,i,itree, tsno(st)%tree(i)%my_tree_number
          stop
       end if
       if(tsno(st)%tree(i)%level.le.0) then
          write(errunit,*) '> Error in check_numbers for tree:',itree,data_tree(itree)%my_number
          write(errunit,*) '> st,i,level:',st,i, tsno(st)%tree(i)%level
          stop
       end if
       
       nstruct_ch = nstruct_ch + 1
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
             i          = nprog
             nbranch_ch = nbranch_ch + 1
          end if
       end if
    end do
    if(nstruct_ch.ne.data_tree(itree)%nb_structs.or.nbranch_ch.ne.data_tree(itree)%nb_branchs) then
       write(errunit,*) '> Error in check_numbers for tree:',itree,data_tree(itree)%my_number
       write(errunit,*) '> nstruct_ch,data_tree(itree)%nb_structs:',nstruct_ch,data_tree(itree)%nb_structs
       write(errunit,*) '> nbranch_ch,data_tree(itree)%nb_branchs:',nbranch_ch,data_tree(itree)%nb_branchs
       stop
    end if

    return

  end subroutine check_numbers

  !*****************************************************************************************************************
  subroutine deallocate_all
    
    implicit none
    integer(kind=4) :: st

    do st = 1, nsteps
       if(tsno(st)%nb_of_structs.gt.0) then
          deallocate(tsno(st)%tree)
       end if
    end do
    deallocate(tsno)
    deallocate(data_tree)
    
    return
    
  end subroutine deallocate_all

 !*****************************************************************************************************************
  subroutine deallocate_tree_only
    
    implicit none
    integer(kind=4) :: st,i
    
    do st = 1, nsteps
       if(tsno(st)%nb_of_structs.gt.0) then
          do i = 1,tsno(st)%nb_of_structs
             if(tsno(st)%tree(i)%ndads.gt.0) then
                deallocate(tsno(st)%tree(i)%dads)
                deallocate(tsno(st)%tree(i)%mass_dads)
             end if
          end do
          deallocate(tsno(st)%tree)
       end if
    end do
    deallocate(tsno)
    
    return
    
  end subroutine deallocate_tree_only
  
  !*****************************************************************************************************************
  
end module init_tree
