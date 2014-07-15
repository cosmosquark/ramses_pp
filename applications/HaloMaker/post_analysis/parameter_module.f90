module parameter_module
  
  ! simple module which holds declarations of run parameters and routine to read them from parameter file

  public

  character(512)  :: repository_jeje ! where output directories of ramses are. 
  integer(kind=4) :: output_nb       ! snapshot number (will read repository_jeje/output_output_nb)
   
  real(kind=8)    :: dnh_fix         ! if >= 0, use as density everywhere (in at/cc)
  real(kind=8)    :: temp_fix        ! if >= 0, use as temperature everywhere (in K)
  real(kind=8)    :: b_fix           ! if > 0, use as Doppler parameter everywhere (km/s)
  real(kind=8)    :: dust2gas_fix    ! if > 0, use as dust to gas (HI) ratio everywhere.
  logical(kind=4) :: static          ! if true, set all velocities to 0

  character(10)   :: region_shape    ! 'cube' or 'sphere' 
  real(kind=8)    :: region_xmin     ! min x in box units, for 'cube' case
  real(kind=8)    :: region_xmax     ! max x in box units, for 'cube' case
  real(kind=8)    :: region_ymin     ! min y in box units, for 'cube' case
  real(kind=8)    :: region_ymax     ! max y in box units, for 'cube' case
  real(kind=8)    :: region_zmin     ! min z in box units, for 'cube' case
  real(kind=8)    :: region_zmax     ! max z in box units, for 'cube' case
  real(kind=8)    :: region_xc       ! x-coord of sphere center, in box units, for 'sphere' case
  real(kind=8)    :: region_yc       ! y-coord of sphere center, in box units, for 'sphere' case
  real(kind=8)    :: region_zc       ! z-coord of sphere center, in box units, for 'sphere' case
  real(kind=8)    :: region_radius   ! sphere radius, in box units, for 'sphere' case
  real(kind=8)    :: region_radius2  ! square of region radius
  real(kind=8)    :: region_size     ! max extent of region (cube or sphere)

  logical(kind=4) :: continuum       ! flat continuum if true, monochromatic if false
  integer(kind=4) :: nphot_tot       ! per frequency (bin)
  real(kind=8)    :: param_v_source  ! emission frequency in km/s
  real(kind=8)    :: v_source_min    ! min emission frequency, in km/s
  real(kind=8)    :: v_source_max    ! max emission frequency, in km/s
  real(kind=8)    :: v_source_step   ! sampling step, in km/s
  character(10)   :: sources         ! 'point' or 'stars'
  real(kind=8)    :: max_stellar_age ! maximum age of emitting stars, in Myr
  real(kind=8)    :: source_xx       ! x-coord of source, in box units
  real(kind=8)    :: source_yy       ! y-coord of source, in box units
  real(kind=8)    :: source_zz       ! z-coord of source, in box units
  
  logical(kind=4) :: recoil          ! set to false to neglect recoil effect
  logical(kind=4) :: deuterium       ! set to false to ignore Deuterium contribution
  integer(kind=4) :: dipol           ! (1 == dipole, 2 == isotropic, 3 == QM, 4 == QM for Lya + Henyey-Greenstein for dust)
  
  integer(kind=4) :: iran            ! seed for random number generation 
  logical(kind=4) :: verbose         ! set to true to print awesome messages at runtime
  
contains
  
  subroutine read_parameter_file
    
    implicit none 

    character(512)  :: line,name,value
    integer(kind=4) :: i

    ! define default values for all parameters 
    repository_jeje       = '/data2/blaizot/Slab/'  ! where output directories of ramses are. 
    output_nb        = 1         ! snapshot number (will read repository/output_output_nb)
    dnh_fix          = 0.0       ! if > 0, use as density everywhere (in at/cc)
    temp_fix         = 0.0       ! if > 0, use as temperature everywhere (in K)
    b_fix            = 0.0       ! if > 0, use as Doppler parameter everywhere (km/s)
    dust2gas_fix     = 0.0       ! if > 0, use as dust to gas (HI) ratio everywhere.
    static           = .false.   ! if true, set all velocities to 0
    region_shape     = 'cube'    ! 'cube' or 'sphere' 
    region_xmin      = 0.0       ! min x in box units, for 'cube' case
    region_xmax      = 1.0       ! max x in box units, for 'cube' case
    region_ymin      = 0.0       ! min y in box units, for 'cube' case
    region_ymax      = 1.0       ! max y in box units, for 'cube' case
    region_zmin      = 0.0       ! min z in box units, for 'cube' case
    region_zmax      = 1.0       ! max z in box units, for 'cube' case
    region_xc        = 0.5       ! x-coord of sphere center, in box units, for 'sphere' case
    region_yc        = 0.5       ! y-coord of sphere center, in box units, for 'sphere' case
    region_zc        = 0.5       ! z-coord of sphere center, in box units, for 'sphere' case
    region_radius    = 0.2       ! sphere radius, in box units, for 'sphere' case
    continuum        = .false.   ! flat continuum if true, monochromatic if false
    nphot_tot        = 100       ! per frequency (bin)
    param_v_source   = 0.0       ! emission frequency in km/s
    v_source_min     = -5000.0   ! min emission frequency, in km/s
    v_source_max     = 5000.0    ! max emission frequency, in km/s
    v_source_step    = 10.0      ! sampling step, in km/s
    sources          = 'stars'   ! 'point' or 'stars'
    max_stellar_age  = 40.0      ! maximum age of emitting stars, in Myr
    source_xx        = 0.5       ! x-coord of source, in box units
    source_yy        = 0.5       ! y-coord of source, in box units
    source_zz        = 0.5       ! z-coord of source, in box units
    recoil           = .true.    ! set to false to neglect recoil effect
    deuterium        = .true.    ! set to false to ignore Deuterium contribution
    dipol            = 4         ! (1 == dipole, 2 == isotropic, 3 == QM, 4 == QM for Lya + Henyey-Greenstein for dust)
    iran             = -100      ! seed for random number generation 
    verbose          = .false.   ! set to true to print awesome messages at runtime

    
    ! read parameters from params.dat file
    open(unit=10,file='params.dat',status='old',form='formatted')
    do
       read (10,'(a)',end=2) line
       i = scan(line,'=')
       if (i==0 .or. line(1:1)=='#') cycle  ! skip blak or commented lines
       name=trim(adjustl(line(:i-1)))
       value=trim(adjustl(line(i+1:)))
       i = scan(value,'!')
       if (i /= 0) value = trim(adjustl(value(:i-1)))
       write(*,'(a,a15,a3,a)') '> ',trim(name),' : ',trim(value)
       select case (trim(name))
       case ('repository_jeje')
          repository_jeje = trim(value)
       case ('output_nb')
          read(value,*) output_nb
       case ('dnh_fix')
          read(value,*) dnh_fix
       case ('temp_fix')
          read(value,*) temp_fix
       case ('b_fix') 
          read(value,*) b_fix
       case ('dust2gas_fix')
          read(value,*) dust2gas_fix
       case ('static')
          read(value,*) static
       case ('shape') 
          region_shape = trim(value)
       case ('xmin') 
          read(value,*) region_xmin 
       case ('xmax') 
          read(value,*) region_xmax
       case ('ymin') 
          read(value,*) region_ymin 
       case ('ymax') 
          read(value,*) region_ymax
       case ('zmin') 
          read(value,*) region_zmin 
       case ('zmax') 
          read(value,*) region_zmax
       case ('xc') 
          read(value,*) region_xc
       case ('yc') 
          read(value,*) region_yc
       case ('zc') 
          read(value,*) region_zc
       case ('radius')
          read(value,*) region_radius
       case ('continuum')
          read(value,*) continuum 
       case ('nphotons')
          read(value,*) nphot_tot
       case ('v_source') 
          read(value,*) param_v_source
       case ('v_source_min')
          read(value,*) v_source_min
       case ('v_source_max')
          read(value,*) v_source_max
       case ('v_source_step') 
          read(value,*) v_source_step
       case ('sources')
          sources = trim(value)
       case ('max_stellar_age')
          read(value,*) max_stellar_age
       case ('x_source','source_xx','xx_em')
          read(value,*) source_xx
       case ('y_source','source_yy','yy_em')
          read(value,*) source_yy
       case ('z_source','source_zz','zz_em')
          read(value,*) source_zz
       case ('recoil')
          read(value,*) recoil
       case ('deuterium')
          read(value,*) deuterium
       case ('dipole_model')
          read(value,*) dipol
       case ('iran')
          read(value,*) iran
       case ('verbose')
          read(value,*) verbose
       end select
    end do
2   close(10)

    ! define some stuff to save computations later on
    region_radius2 = region_radius * region_radius
    select case (trim(region_shape))
    case ('cube')       
       region_size = max(region_xmax - region_xmin,region_ymax - region_ymin,region_zmax - region_zmin)
    case ('sphere')
       region_size = region_radius * 2.0d0
       ! following defs are useful for first pass with hilbert curve ... (second pass used function in_region below). 
       region_xmin = region_xc - region_radius
       region_xmax = region_xc + region_radius
       region_ymin = region_yc - region_radius
       region_ymax = region_yc + region_radius
       region_zmin = region_zc - region_radius
       region_zmax = region_zc + region_radius
    case default
       write(*,*) 'ERROR in parameter_module:read_parameter_file : regrion_shape unknown'
       write(*,*) trim(region_shape)
       stop
    end select
    
    return

  end subroutine read_parameter_file


  function in_region(x,y,z)

    ! tests if a point (x,y,z, in box units) is in the selected region
    ! -> returns true if so, false otherwise
    
    implicit none

    logical(kind=4)         :: in_region
    real(kind=8),intent(in) :: x,y,z
    real(kind=8)            :: d2

    in_region = .true. 
    select case (trim(region_shape))
    case ('cube') 
       if ( x.lt.region_xmin .or. x.gt.region_xmax .or. & 
            y.lt.region_ymin .or. y.gt.region_ymax .or. &
            z.lt.region_zmin .or. z.gt.region_zmax) then
          in_region = .false.
       end if
    case ('sphere')
       d2 = (x - region_xc)**2 + (y - region_yc)**2 + (z - region_zc)**2
       if (d2 > region_radius2) in_region = .false.
    case default
       write(*,*) 'ERROR in parameter_module:in_region : regrion_shape unknown'
       stop
    end select

    return

  end function in_region

end module parameter_module









