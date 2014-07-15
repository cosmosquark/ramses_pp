module tree_defs

public

  type tree_info
     integer(kind=4)          :: my_number        ! id of halo in the simualtion
     integer(kind=4)          :: my_timestep      ! step of halo in the simulation
     integer(kind=4)          :: my_tree_number   ! id of last halo progenitor
     real(kind=4)             :: m                ! halo mass 1e11 M_sun
     real(kind=4)             :: mvir             ! halo virial mass 1e11 M_sun
     integer(kind=4)          :: level            ! level in the halo structure tree
     integer(kind=4)          :: hosthalo         ! id of the host halo
     integer(kind=4)          :: hostsub          ! id of the direct host (either halo or subhalo)
     integer(kind=4)          :: nextsub          ! id of the next subhalo in the halo structure tree 
                                                  ! (linked list of all subhalos for a given halo starting from %hosthalo)
     integer(kind=4)          :: ndads            ! nb of progenitors
     integer(kind=4),pointer  :: dads(:)          ! id of progentors 
     real(kind=4),pointer     :: mass_dads(:)     ! mass fraction obtained from progenitors (M(prog to halo)/M(halo))
     integer(kind=4)          :: son              ! main descendant
     integer(kind=4)          :: fst_prog         ! main progenitor
     integer(kind=4)          :: nxt_prog         ! next secondary progenitor 
                                                  ! (linked list of all progenitors of a given halo, strating from %fst_prog)
     integer(kind=4)          :: frag             ! fragment flag, 0 if normal, 1 if fragment, 2 if all progenitors are fragments
  end type tree_info

  type time_info
     real(kind=4)             :: aexp             ! expansion factor
     real(kind=4)             :: age_univ         ! age of univers Gyr
     real(kind=4)             :: redshift         ! redshift 
     real(kind=4)             :: omega_t          ! omega matter
     integer(kind=4)          :: nb_of_halos      ! nb of halos in the tree file for this step
     integer(kind=4)          :: nb_of_subhalos   ! nb of subhalos in the tree file for this step
     integer(kind=4)          :: nb_of_structs    ! nb of halos and subhalos in the tree file for this step
     integer(kind=4), pointer :: halo_number(:)   ! halo number in the whole simulation
     type(tree_info), pointer :: tree(:)          ! merger trees
  end type time_info

  type tree_data
     integer(kind=4)          :: my_number        ! corresponding final halo number in the simulation 
     integer(kind=4)          :: nb_structs       ! nb of halos and subhalos in the merger tree
     integer(kind=4)          :: nb_branchs       ! nb of branches in the merger tree
     integer(kind=4)          :: nb_frags         ! nb of fragments in the merger tree
     integer(kind=4)          :: nb_target        
  end type tree_data

  type tree_out_info
     integer(kind=4)          :: my_halo_number   ! id of halo in the simulation
     integer(kind=4)          :: my_step_number   ! st of halo in the simulation
     integer(kind=4)          :: my_tree_number   ! id of the merger tree in the whole simulation
     integer(kind=4)          :: my_branch_number ! id of the branch in the full merger tree
     integer(kind=4)          :: my_tree_id       ! tree index in the tree_file
     integer(kind=4)          :: my_halo_id       ! halo index in the tree_file
     integer(kind=4)          :: my_branch_id     ! branch index in the zoomed merger tree
     integer(kind=4)          :: stend            ! last step for tree zoomed visualisation
     integer(kind=4)          :: ststart          ! first step for tree zoomed visualisation
     integer(kind=4)          :: flag_out         ! 0 if not found, 1 if found in file, -1 if already outputed
  end type tree_out_info

  type(time_info),allocatable     :: tsno(:)
  type(tree_data),allocatable     :: data_tree(:)
  type(tree_out_info),allocatable :: tree_list(:)
  
  integer(kind=4)    :: nsteps
  character(len=3)   :: method
  character(len=200) :: file_tree
  character(len=200) :: dir_tree
  character(len=200) :: dir_out
  integer(kind=4)    :: n_tree_files
  integer(kind=4)    :: n_tree_out,n_tree_left,n_tree_file_out
  integer(kind=4)    :: nbranchmax 
  real(kind=4)       :: res_force
  
  !======================================================================
  ! Definitions specific to input/output
  !======================================================================
  character(80)      :: data_dir
  integer(kind=4)    :: errunit = 0
  !======================================================================
  
end module tree_defs
