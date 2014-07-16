module jeje 
  
  use parameter_module

  public

#ifdef JEJE
  
  ! define model for polytrope (see end of ramses2real)
  ! 1. (polytrope_model == 1) -> temperature is set to 1.d4K and gas to fully neutral (nHI = nHtot)
  ! 2. (polytrope_model == 2) -> density is set to zero in high & cold regions so they become transparent
  ! 3. (polytrope_model == 3) -> polytrope for GALAXY_CLUMP_WINDS (T<100K, n>10 at/cm3)
  ! 4. (polytrope_model == 4) -> temperature is set to 1.d4K and gas to fully neutral (nHI = nHtot) for nhtot>0.3
  !  integer(kind=4),parameter   :: polytrope_model = 3
  integer(kind=4),parameter   :: polytrope_model = 4
  ! define a model for self-shielding ... 
  ! 1. (selfshielding_model = 0) -> no model (take neutral fraction from RAMSES, assuming PIE)
  ! 2. (selfshielding_model = 1) -> assume gas to be neutral above nh = 0.1 at / cc
  integer(kind=4),parameter   :: selfshielding_model = 0

  type cooling_table
     integer(kind=4)          :: n11
     integer(kind=4)          :: n22
     real(kind=8),allocatable :: nH(:)    
     real(kind=8),allocatable :: T2(:)    
     real(kind=8),allocatable :: T2eq(:)    
     real(kind=8),allocatable :: metal(:,:)  
     real(kind=8),allocatable :: cool(:,:)  
     real(kind=8),allocatable :: heat(:,:)  
     real(kind=8),allocatable :: metal_prime(:,:)  
     real(kind=8),allocatable :: cool_prime(:,:)  
     real(kind=8),allocatable :: heat_prime(:,:)  
     real(kind=8),allocatable :: mu(:,:)  
     real(kind=8),allocatable :: spec(:,:,:)  ! last dimension (6) is n_e, n_HI, n_HII, n_HeI, n_HeII, n_HeIII
  end type cooling_table
  type(cooling_table) :: cooling

  type cool_interp
     integer(kind=4)  :: n_nH
     real(kind=8)     :: nH_start,nH_step
     integer(kind=4)  :: n_T2
     real(kind=8)     :: T2_start,T2_step
  end type cool_interp
  type(cool_interp)  :: cool_int

  ! public conversion factors -> see subroutine get_conversion_factors
  real(kind=8)   :: dp_scale_l 
  real(kind=8)   :: dp_scale_d 
  real(kind=8)   :: dp_scale_t 
  real(kind=8)   :: dp_scale_nH   
  real(kind=8)   :: dp_scale_v    
  real(kind=8)   :: dp_scale_T2   
  real(kind=8)   :: dp_scale_zsun 
  
  real(kind=8),parameter :: XH      = 0.76                ! mass fraction of H 
  real(kind=8),parameter :: mH      = 1.6600000d-24
  real(kind=8),parameter :: kB      = 1.3806200d-16

  ! USAGE 
  !
  ! 1/ call init_cooling (once per output number) 
  ! 2/ call ramses2real 
  ! 

contains

!*****************************************************************************************************************

  subroutine init_cooling

    implicit none 
    
    call get_conversion_factors
#ifndef RT
    call read_cooling
#endif

  end subroutine init_cooling

!*****************************************************************************************************************
  
  subroutine get_conversion_factors

    implicit none

    character(512)             :: filename
    integer(kind=4)            :: i
    
    ! read scales from info_xxx.txt file
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository_jeje),'/output_',output_nb,"/info_",output_nb,".txt"        
    open(unit=44,file=filename,form='formatted',status='old')
    do i = 1, 15    ! skip useless header
       read(44,*)
    end do
    read(44,'(13X,E23.15)')dp_scale_l
    read(44,'(13X,E23.15)')dp_scale_d
    read(44,'(13X,E23.15)')dp_scale_t
    close(44)
    dp_scale_nH   = XH/mH * dp_scale_d      ! convert mass density (code units) to numerical density of H atoms
    dp_scale_v    = dp_scale_l/dp_scale_t   ! converts velocities (code units) to cm/s
    dp_scale_T2   = mH/kB * dp_scale_v**2   ! converts T(/mu) into K
    dp_scale_zsun = 1.d0/0.0127             ! solar metallicity

    return

  end subroutine get_conversion_factors

!*****************************************************************************************************************

  subroutine read_cooling

    implicit none

    character(1024)            :: filename
    integer(kind=4)            :: n1,n2

    ! initialize cooling variables
    call clear_cooling

    ! read cooling variables from cooling.out file
    write(filename,'(a,a,i5.5,a,i5.5,a)') trim(repository_jeje),'/output_',output_nb,"/cooling_",output_nb,".out"
    open(unit=44,file=filename,form='unformatted')
    read(44) n1,n2
    cooling%n11 = n1
    cooling%n22 = n2
    allocate(cooling%nH(n1),cooling%T2(n2),cooling%T2eq(n1))
    allocate(cooling%cool(n1,n2),cooling%heat(n1,n2),cooling%mu(n1,n2))
    allocate(cooling%cool_prime(n1,n2),cooling%heat_prime(n1,n2),cooling%metal_prime(n1,n2))
    allocate(cooling%metal(n1,n2),cooling%spec(n1,n2,6))
    read(44)cooling%nH
    read(44)cooling%T2
    read(44)cooling%T2eq
    read(44)cooling%cool
    read(44)cooling%heat
    read(44)cooling%metal
    read(44)cooling%cool_prime
    read(44)cooling%heat_prime
    read(44)cooling%metal_prime
    read(44)cooling%mu
    read(44)cooling%spec
    close(44)
    
    ! define useful quantities for interpolation 
    cool_int%n_nh     = n1
    cool_int%nh_start = minval(cooling%nh)
    cool_int%nh_step  = cooling%nh(2) - cooling%nh(1)
    cool_int%n_t2     = n2
    cool_int%t2_start = minval(cooling%t2)
    cool_int%t2_step  = cooling%t2(2) - cooling%t2(1)
    
    return

  end subroutine read_cooling

!*****************************************************************************************************************

  subroutine clear_cooling

    implicit none

    if (cooling%n11 > 0 .or. cooling%n22 > 0) then 
       cooling%n11 = 0
       cooling%n22 = 0
       deallocate(cooling%nH,cooling%T2,cooling%T2eq,cooling%metal)
       deallocate(cooling%cool,cooling%heat,cooling%metal_prime,cooling%cool_prime)
       deallocate(cooling%heat_prime,cooling%mu,cooling%spec)
    end if
    
    cool_int%n_nh = 0
    cool_int%nh_start = 0.0d0
    cool_int%nh_step  = 0.0d0
    cool_int%n_t2 = 0
    cool_int%t2_start = 0.0d0
    cool_int%t2_step  = 0.0d0
    
    return

  end subroutine clear_cooling

!*****************************************************************************************************************
#ifdef RT
  subroutine ramses2real(ncell,temp,nh,vx,vy,vz,ionfrac)
#else
  subroutine ramses2real(ncell,temp,nh,vx,vy,vz)
#endif
    
    implicit none 

    ! converts code units to : 
    ! - velocity in cm/s
    ! - neutral H density in at / cm^3
    ! - temperature (no mu!) in K
    ! -> inputs (temp, nh, etc) are assumed to be as read from the snapshot (except temp=var(5)/var(1))

    integer(kind=4),intent(in) :: ncell
    real(kind=8),intent(inout) :: nh(ncell),temp(ncell),vx(ncell),vy(ncell),vz(ncell)
#ifdef RT
    real(kind=8),intent(in)    :: ionfrac(3,ncell)
    real(kind=8)               :: mu
#endif
    integer(kind=4) :: ihx,ihy,i
    real(kind=8)    :: xx,yy,dxx1,dxx2,dyy1,dyy2,f,nhtot
    integer(kind=4) :: if1,if2,jf1,jf2

    vx   = vx * dp_scale_v ! cm/s
    vy   = vy * dp_scale_v ! cm/s
    vz   = vz * dp_scale_v ! cm/s
    
    temp = (temp/nh) * dp_scale_T2  ! convert P/rho (code units) into T/mu (in Kelvin)
    nh   = nh * dp_scale_nH    ! number density (/cc) of hydrogen atoms (neutral or not)

#ifdef RT

    ! if RT, compute mu directly from ionized fractions of H and He and derive T
    ! compute neutral H directly from ionized H fraction. 
    do i = 1,ncell
       if (nh(i) > 0.0d0) then  ! only deal with leaf cells
          nhtot   = nh(i)
          mu      = 1.d0/( XH*(1.d0+ionfrac(1,i)) + 0.25*(1.d0-XH)*(1.d0+ionfrac(2,i)+2.*ionfrac(3,i)) )  ! assumes no metals.
          temp(i) = temp(i) * mu
          nh(i)   = nh(i) * ionfrac(1,i)
          ! Deal with the polytrope
          select case (polytrope_model)
          case (1)
             ! simply force the temperature back down to 1d4 K and gas to be neutral
             if (nhtot > 0.1 .and. temp(i) > 1.d4) then 
                temp(i) = 1.d4
                nh(i)   = nhtot
             end if
          case (2)
             ! make polytrope cells (ie. ISM) transparent
             if (nhtot > 0.1) then 
                nh(i)   = 0.0d0
             end if
          case (3)
             if (nhtot > 10.0d0 .and. temp(i) > 1.d2) then 
                temp(i) = 1.d2
                nh(i)   = nhtot
             end if
          case default
             write(6,*) 'polytrope model not known'
          end select
       end if

       end if
    end do

#else

    ! if no RT, use tabulated equilibrium solutions to infer mu, T, and nHI
    do i = 1,ncell 

       if (nh(i) > 0.0d0) then  ! check whether density is > 0 as a proxy for testing whether cell is a leaf... 
          if (temp(i) <=0) then 
             !print*,'Oh nooo!',temp(i),nh(i)
             !stop
             cycle
          end if

          xx  = log10(nh(i))
          ihx = int((xx - cool_int%nh_start)/cool_int%nh_step) + 1
          if (ihx < 1) then 
             ihx = 1 
          else if (ihx > cool_int%n_nh) then
             ihx = cool_int%n_nh
          end if
          yy  = log10(temp(i))
          ihy = int((yy - cool_int%t2_start)/cool_int%t2_step) + 1
          if (ihy < 1) then 
             ihy = 1 
          else if (ihy > cool_int%n_t2) then
             ihy = cool_int%n_t2
          end if
          ! 2D linear interpolation:
          if (ihx < cool_int%n_nh) then 
             dxx1  = max(xx - cooling%nh(ihx),0.0d0) / cool_int%nh_step 
             dxx2  = min(cooling%nh(ihx+1) - xx,cool_int%nh_step) / cool_int%nh_step
             if1  = ihx
             if2  = ihx+1
          else
             dxx1  = 0.0d0
             dxx2  = 1.0d0
             if1  = ihx
             if2  = ihx
          end if
          if (ihy < cool_int%n_t2) then 
             dyy1  = max(yy - cooling%t2(ihy),0.0d0) / cool_int%t2_step
             dyy2  = min(cooling%t2(ihy+1) - yy,cool_int%t2_step) / cool_int%t2_step
             jf1  = ihy
             jf2  = ihy + 1
          else
             dyy1  = 0.0d0
             dyy2  = 1.0d0
             jf1  = ihy
             jf2  = ihy
          end if
          if (abs(dxx1+dxx2-1.0d0) > 1.0d-6 .or. abs(dyy1+dyy2-1.0d0) > 1.0d-6) then 
             write(*,*) 'Fucked up the interpolation ... '
             print*,dxx1+dxx2,dyy1+dyy2
             stop
          end if
          
          ! neutral H density 
          f = dxx1 * dyy1 * cooling%spec(if2,jf2,2) + dxx2 * dyy1 * cooling%spec(if1,jf2,2) &
               & + dxx1 * dyy2 * cooling%spec(if2,jf1,2) + dxx2 * dyy2 * cooling%spec(if1,jf1,2)
          nhtot = nh(i)
          nh(i) = real(10.0d0**f,4)   ! nHI (cm^-3)
          ! GET MU to convert T/MU into T ... 
          f = dxx1 * dyy1 * cooling%mu(if2,jf2) + dxx2 * dyy1 * cooling%mu(if1,jf2) &
               & + dxx1 * dyy2 * cooling%mu(if2,jf1) + dxx2 * dyy2 * cooling%mu(if1,jf1)
          temp(i) = temp(i) * f   !! This is now T (in K) with no bloody mu ... 

          ! Deal with self-shielded regions
          if (selfshielding_model == 1) then
             ! -> assume that gas at densities > 1.d-2 at/cc and T < 3d4 is self-shielded and hence fully neutral. 
             if (nhtot > 0.1 .and. temp(i) < 3.d4) then 
                nh(i) = nhtot
             end if
          end if
          
          ! Deal with the polytrope in dense regions
          select case (polytrope_model)
          case (1)
             ! simply force the temperature back down to 1d4 K 
             if (nhtot > 0.1 .and. temp(i) > 1.d4) then 
                temp(i) = 1.d4
                nh(i)   = nhtot
             end if
          case (2)
             ! make polytrope cells (ie. ISM) transparent)
             if (nhtot > 0.1 .and. temp(i) > 1.d4) then 
                nh(i)   = 0.0d0  
             end if
          case (3)
             if (nhtot > 10.0 .and. temp(i) > 1.d2) then 
                nh(i)   = nhtot 
                temp(i) = 1.d2
             end if
          case (4)
             ! simply force the temperature back down to 1d4 K 
             if (nhtot > 0.3 .and. temp(i) > 1.d4) then 
                temp(i) = 1.d4
                nh(i)   = nhtot
             end if
          case default
             write(6,*) 'polytrope model not known'
          end select
       end if
    end do
    
#endif

 
  end subroutine ramses2real

!*****************************************************************************************************************

#endif

end module jeje




