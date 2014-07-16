module ramses_util
implicit none

contains

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif

end subroutine title

!================================================================
Subroutine read_info(arg,scale,scale_l,scale_d,scale_t,t,aexp,omega_m,omega_l,omega_k)
  implicit none
  integer::ncpu,ndim,npart,ilun
  integer::nlevelmax,levelmin
  integer::ngridmax,nstep_coarse
  real(KIND=8)::t,aexp,scale,h0
  real(KIND=8)::omega_m,omega_l,omega_k,omega_b
  real(KIND=8)::scale_l,scale_d,scale_t
  logical::ok
  character*5::nchar
!  character*128::nomfich,nomdir
  character*256::nomfich,nomdir
  integer::iargc,ipos
  character(len=128)::arg


  nomdir=trim(arg)//'/'
  ipos=INDEX(nomdir,'output_')
  nchar=nomdir(ipos+7:ipos+13)
  ! Lecture du fichier info   
                                                                                                      
  nomfich=TRIM(nomdir)//'info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok)
  if (.not. ok) then
     write(*,*)'File '//TRIM(nomfich)//' not found'
  else
     ilun=10
     open(unit=ilun,file=nomfich,status='old',form='formatted')
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
     close(ilun)
  endif

end subroutine read_info

!================================================================
Subroutine read_info_nstep(arg,nstep_coarse,aexp,scale_l,scale_d,scale_t)
  implicit none
  integer::ncpu,ndim,npart,ilun
  integer::nlevelmax,levelmin
  integer::ngridmax,nstep_coarse
  real(KIND=8)::t,aexp,scale,h0
  real(KIND=8)::omega_m,omega_l,omega_k,omega_b
  real(KIND=8)::scale_l,scale_d,scale_t
  logical::ok
  character*5::nchar
  character*50::nomfich,nomdir
  integer::iargc,ipos
  character(len=128)::arg


  nomdir=trim(arg)//'/'
  ipos=INDEX(nomdir,'output_')
  nchar=nomdir(ipos+7:ipos+13)
  ! Lecture du fichier info   
                                                                                                      
  nomfich=TRIM(nomdir)//'info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok)
  if (.not. ok) then
     write(*,*)'File '//TRIM(nomfich)//' not found'
  else
     ilun=10
     open(unit=ilun,file=nomfich,status='old',form='formatted')
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
     close(ilun)
  endif

end subroutine read_info_nstep

!================================================================
!================================================================
subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min,age_tot
  real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
!  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

!!$  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
!!$     write(*,*)'Error: non-physical cosmological constants'
!!$     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
!!$     write(*,*)'The sum must be equal to 1.0, but '
!!$     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
!!$     stop
!!$  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0
  
  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau
     
     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
  end do

  age_tot=-t
  write(*,666)-t
  666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable
  
  axp_t = 1.d0
  t = 0.d0
  axp_tau = 1.d0
  tau = 0.d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if

  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman

function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau

function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt
!================================================================
!================================================================
 
end module ramses_util


