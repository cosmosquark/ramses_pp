module utils

  ! include general routines which have to be moved here because used by subbox and compute_halo_props

  use halo_defs

  public

contains

!***********************************************************************

  subroutine correct_for_periodicity(dr)

  ! subroutine corrects for the fact that if you have periodic boundary conditions,
  ! then groups of particles that sit on the edge of the box can have part of their
  ! particles on one side and part on the other. So we have to take out a box-length
  ! when measuring the distances between group members if needed.

    implicit none

    type (vector) :: dr

    if (FlagPeriod == 0) return  !--> NO PERIODIC BCs 
    
    if (dr%x > + Lbox_pt2) dr%x = dr%x - Lbox_pt
    if (dr%x <= - Lbox_pt2) dr%x = dr%x + Lbox_pt 

    if (dr%y > + Lbox_pt2) dr%y = dr%y - Lbox_pt
    if (dr%y <= - Lbox_pt2) dr%y = dr%y + Lbox_pt 

    if (dr%z > + Lbox_pt2) dr%z = dr%z - Lbox_pt
    if (dr%z <= - Lbox_pt2) dr%z = dr%z + Lbox_pt 

    return

  end subroutine correct_for_periodicity

!***********************************************************************

  subroutine correct_for_periodicity_code_units(x)
    
    ! same as correct_for_periodicity but argument is only one coord in code units (i.e. from -0.5 to 0.5)

    implicit none

    real(kind=4) :: x

    if (x >  0.5d0) then 
       x = x - 1.0d0
    end if
    if (x <= -0.5d0) then 
       x = x + 1.0d0
    end if

    return

  end subroutine correct_for_periodicity_code_units

!***********************************************************************

  subroutine title(n,nchar)
  
    implicit none

    integer(kind=4) :: n
    character*5     :: nchar
    character*1     :: nchar1
    character*2     :: nchar2
    character*3     :: nchar3
    character*4     :: nchar4
    character*5     :: nchar5
    

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
  
!***********************************************************************                

end module utils



