!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! @author David Sullivan (D.Sullivan@sussex.ac.uk)
!
! compile with: mpif90 ramses_part2ascii.f90 -o ramses_part2ascii
! TODO Speed up writing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ramses_part2ascii

	use mpi

	implicit none

	integer				:: i, j, rank, icpu, ncpu, nfiles, mpi_reduce_buffer
	integer 			:: firstfile, lastfile, files_per_cpu, ifile
	integer				:: mpi_ierr
	character(len=5)	:: dir_number, suffix
	real(kind=8)		:: x, y, z, vx, vy, vz, m, u, m_group
	integer (kind=4)	:: npart, ngaspart, nramsespart, ndmpart, nstarpart, nsink, intdummy
	integer				:: part_file, info_file, output_file

	!arrays
	real(kind=8),dimension(:,:), allocatable :: ramsespart_pos, dmpart_pos, starpart_pos
	real(kind=8),dimension(:,:), allocatable :: ramsespart_vel, dmpart_vel, starpart_vel
	real(kind=8),dimension(:), allocatable :: ramsespart_m, dmpart_m, starpart_m
	real(kind=8),dimension(:), allocatable :: ramsespart_id, dmpart_id, starpart_id

	! variables for file management
	character(len=128) :: dir_name, out_dir, out_name, part_filename, info_filename

	!bools
	logical	:: verbose = .true.
	logical :: debug = .true.

	!Init MPI environment
	call MPI_INIT(mpi_ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD, icpu, mpi_ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, mpi_ierr)

	call read_args(dir_name, out_dir)
	if (icpu == 0) call hello(dir_name, out_dir, icpu, ncpu)

	!Construct filename, inquire existence and assign I/O unit number
	dir_number=dir_name( index(dir_name,'output_')+7 : index(dir_name,'output_')+13 ) 
	!dir_number = '00011'

	if (icpu == 0) then

		info_filename   = trim(dir_name) // '/info_'   // trim(dir_number) // '.txt'
		call inquire_file (info_filename)
		info_file   = 988

		write(*,*) 'Info filename: ', info_filename

		open(unit=info_file,file=info_filename,form='formatted',status='old',action='read', &
			& position = 'rewind')

		!Read number of files from the info file
		read(info_file,'(13X,I11)')	nfiles

		!At the moment this is restricted to nfiles == ncpu
		if (.not. debug .and. nfiles /= ncpu ) call terminate(1)

		! Broadcast nfiles
		call MPI_Bcast(nfiles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_ierr)

	end if

	files_per_cpu = ( nfiles+(ncpu-1) ) / ncpu
  	firstfile = (icpu*files_per_cpu) + 1
  	lastfile  = min ( (icpu+1)*files_per_cpu, nfiles )

	!Write the output file number based on this MPI rank
	write(suffix,'(i5.5)') (icpu+1)
	part_filename =   trim(dir_name) // '/part_'  // trim(dir_number) // '.out' // trim(suffix)
	!Assign a unique I/O unit number
	part_file    = 10*icpu + 1

	!Inquire if the particle file exists
	call inquire_file(part_filename)

	write (*,*) 'CPU ', icpu, 'working on file ', part_filename

	!Read the particle file
	open(unit=part_file,file=part_filename,status='old',form='unformatted',action='read')

	call read_particle_file

	!We're done, clode the stream
	close(part_file)

	!Assign the same I/O units number across all MPI ranks as we want to write to a single file (as this was designed for running powmes)
	output_file = 990
	out_name = trim(out_dir) // '/asciitest.ascii'
	!open(unit=output_file,file=out_name,action="write")

	if (icpu == 0) then ! Write number of particles
		write(*,*) 'Total particles: ', npart
		open(unit=output_file,file=out_name,action="write",status="replace")
		write(output_file,*) npart
		close(output_file)
	end if

	!I do not like this...
	!Have eack rank append onto the file one at a time
	do rank = 0, ncpu
		if (icpu == rank) then
			write(*,*) 'CPU on rank ', icpu, ' writing positions and mass'
			open(unit=output_file,file=out_name,action="write",status="old",position='append')
			do i = 1, nramsespart
				write(output_file,*) ramsespart_pos(0, i), ramsespart_pos(1, i), ramsespart_pos(2, i), ramsespart_m(i)
			end do
			close(output_file)
		end if
		call MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
	end do

	!Close up
	close(output_file)
	call MPI_FINALIZE(mpi_ierr)

	! Dealocate arrays and wrap up!
	deallocate(ramsespart_pos)
	deallocate(ramsespart_vel)
	deallocate(ramsespart_m)

	if (icpu == 0) write(*,*) 'Done'

	contains

	!---------------------------------------------------------------------------------------!

	subroutine read_particle_file

		read(part_file) ! skip nfiles
		read(part_file) ! skip ndim
		read(part_file)	nramsespart
		read(part_file) ! skip local seed
		read(part_file) nstarpart
		read(part_file) ! skip mstar_tot
		read(part_file) ! skip mstar_lost
		read(part_file) nsink

		!Sum all particle counts
		call MPI_ALLREDUCE(nramsespart,mpi_reduce_buffer, &
	       & 1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)

		npart = mpi_reduce_buffer

		!We don't want sink particles...
		if (nsink /= 0) then
			if (icpu == 0) write(*,*) 'Error: Sink particles detected. Aborting'
			call terminate(2)
		end if

		! allocate data arrays for dm and star (=ramses) particles
		allocate ( ramsespart_pos(1:3,1:nramsespart) )
		allocate ( ramsespart_vel(1:3,1:nramsespart) )
		allocate ( ramsespart_m(1:nramsespart) )
		allocate ( ramsespart_id(1:nramsespart) )

		! Read the positions
		do i = 1, 3 ! Loop over dimensions
			read(part_file) (ramsespart_pos(i,j), j=1,nramsespart)
		end do
		!Read the velocities
		do i = 1, 3
			read(part_file) (ramsespart_vel(i,j), j=1,nramsespart)
		end do
		!Read the masses
		read(part_file) ramsespart_m
		!read(part_file) ramsespart_id

	end subroutine read_particle_file

end program ramses_part2ascii

subroutine read_args(dir_name, out_dir)

	implicit none

	character(len=128)				:: arg1
	character(len=128)				:: arg2
	character(len=128),intent(out)	:: dir_name, out_dir

	call getarg(1,arg1)
	call getarg(2,arg2)

	dir_name=trim(arg1)
	out_dir=trim(arg2)

	return

end subroutine read_args

subroutine hello(dir_name, out_dir, icpu, ncpu)

	integer, intent(in)				:: icpu, ncpu
	character(len=128), intent(in)	:: dir_name, out_dir

	write(*,*) 'Rank ', icpu, ' of ', ncpu
	write(*,*) 'Working on dir ', trim(dir_name), ' and writing to ', trim(out_dir)

end subroutine hello

subroutine inquire_file(filename)

	! this checks if a file is present

	implicit none

	character(len=128), intent(in)::filename
	logical :: ok

	inquire(file=filename, exist=ok) 
	if ( .not. ok ) then
		write (*,*) 'Error in input files: ' // trim(filename)// ' not found.'
		call terminate(6)
	end if

	return
end subroutine inquire_file

subroutine terminate(ierr)

	integer,intent(in)	:: ierr
	integer				:: mpi_ierr

	write(*,*) 'Received error signal: ', ierr, 'terminating'
	stop

end subroutine terminate
