MODULE ramses_info
  !--------------------------------------------------------------------------
  ! Module to read and contain the information contained in 
  ! info_xxxxx.txt and info_rt_xxxxx.txt (if it exists).
  ! J. Rosdahl, Jan 2011, modified from code by R. Teyssier.
  !--------------------------------------------------------------------------

  implicit none
  integer::ncpu,ndim,npart
  integer::nlevelmax,levelmin
  integer::ngridmax,nstep_coarse
  real::t,aexp,scale,h0
  real::omega_m,omega_l,omega_k,omega_b
  real::scale_l,scale_d,scale_t
  character(LEN=80)::ordering

  integer::nHvar, nRTvar, nIons, nPacs
  real::X_frac, Y_frac
  real::unit_Np, unit_Fp, rt_c_frac

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE read_info(nsnap)

! Read info files
!------------------------------------------------------------------------
  character(len=128)::filename,dirname
  character(len=128)::nsnap,ifmt,dfmt,varname
  character*5::nchar
  integer::iargc,ipos
  logical::ok
  integer::ilun
  real::value
!------------------------------------------------------------------------

  dirname=trim(nsnap)//'/'
  ipos=INDEX(dirname,'output_')
  nchar=dirname(ipos+7:ipos+13)


  ! Read info file                                                                                                      
  filename=TRIM(dirname)//'info_'//TRIM(nchar)//'.txt'
  inquire(file=filename, exist=ok)
  if (.not. ok) then
     write(*,*)'File '//TRIM(filename)//' not found'
     STOP
  else
     ilun=10
     open(unit=ilun,file=filename,status='old',form='formatted')
     read(ilun,'("ncpu        =",I11)')ncpu
     read(ilun,'("ndim        =",I11)')ndim
     read(ilun,'("levelmin    =",I11)')levelmin
     read(ilun,'("levelmax    =",I11)')nlevelmax
     read(ilun,'("ngridmax    =",I11)')ngridmax
     read(ilun,'("nstep_coarse=",I11)')nstep_coarse
     read(ilun,*)
     read(ilun,'("boxlen      =",E23.15)')scale
     read(ilun,'("time        =",E23.15)')t
     read(ilun,'("aexp        =",E23.15)')aexp
     read(ilun,'("H0          =",E23.15)')h0
     read(ilun,'("omega_m     =",E23.15)')omega_m
     read(ilun,'("omega_l     =",E23.15)')omega_l
     read(ilun,'("omega_k     =",E23.15)')omega_k
     read(ilun,'("omega_b     =",E23.15)')omega_b
     read(ilun,'("unit_l      =",E23.15)')scale_l
     read(ilun,'("unit_d      =",E23.15)')scale_d
     read(ilun,'("unit_t      =",E23.15)')scale_t
     read(ilun,*)
     read(ilun,'("ordering type=",A80)'),ordering
     !read(10,'(A14,A80)')tmp,ordering
     close(ilun)
  endif

#ifdef RT 
  ! Try to read rt info file
  filename=TRIM(dirname)//'info_rt_'//TRIM(nchar)//'.txt'
  inquire(file=filename, exist=ok)
  nHvar   = 6          ; nRTvar  = 0             
  nIons   = 3          ; nPacs   = 0
  X_frac  = 0.76       ; Y_frac  = 0.24
  unit_np = 1.         ; unit_fp = 1.
  rt_c_frac=1.
  if (.not. ok) then
     write(*,*)'File '//TRIM(filename)//' not found'
  else
     open(unit=ilun,file=filename,status='old',form='formatted')
     ifmt='(A13, I30)'
     dfmt='(A13, D0)'
     varname='X_fraction'
     call read_real(ilun, varname, value, 6.)
     nHvar=value
     close(ilun)
  endif
#endif

END SUBROUTINE read_info  

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE read_real(lun, param_name, value, default_value)

! Read a parameter from lun and give it a default value if not found
!------------------------------------------------------------------------
  integer::lun
  character(*)::param_name
  character(128)::line,tmp
  real::value,default_value
!------------------------------------------------------------------------
  rewind(unit=lun)
  line='' ; tmp=''
  do
     read(lun, '(A128)', end=222) line
     if(index(line,trim(param_name)) .eq. 1) then
        print*,line
        print*,'found it '
        read(line,'(A13,E22.7)') tmp, value
        print*,'value=',value
        return
     endif
  end do
  !print,param_name,value
222 value = default_value         ! eof reached, didn't find the parameter

END SUBROUTINE read_real

END MODULE ramses_info
