#
# Makefile for TreeMaker
#
# flags for f90 compilation on crozes
#

# Possible options: 
# -----------------
# - BIG_RUN : for a normal simulation (as opposed to resimulation)

#MACHINE = "meso"
COMPILO = "ifort"

OPTIONS = -DBIG_RUN 

ifeq ($(COMPILO), "pgf90") 
F90C    = pgf90
DEBUG   = -Minform=warn -g -C -Mnofprelaxed
OPTIM   = -O2 -Mnofprelaxed
FFLAGS  = -tp k8-64 -mcmodel=medium -Mlarge_arrays -Mpreprocess -Mextend ${OPTIM}
endif	

ifeq ($(COMPILO), "ifort") 
F90C    = ifort
DEBUG   = -g -CB -fpe0 -traceback -warn all -warn unused #-pg
OPTIM   = -O2 -fp-model precise 
FFLAGS  =  -cpp ${DEBUG}
endif	

ifeq ($(COMPILO), "g95") 
F90C    = g95
DEBUG   = -g -Wall -ftrace=full
OPTIM   = -O2
FFLAGS  = ${DEBUG}
endif	

ifeq ($(COMPILO)$,"gfortran")
# option when compiling with pgf90 (meso, horizon1)
F90C    = gfortran
DEBUG   = -g -Wall -ftrace=full
OPTIM   = -O2
FFLAGS  =  $(OPTIM)
endif

.SUFFIXES: .o .f90 
# Rules:
%.o: %.f90
	$(F90C) $(FFLAGS) $(OPTIONS) -c $*.f90

default:        Show_Trees clean


#
# erase useless files:
#

clean: 	
	rm -f *.o
	rm -f *.mod
	rm -f *~

Show_Trees:  tree_defs.o init_tree.o input_output.o out_tree_mod.o show_trees.o
	${F90C} ${FFLAGS} ${OPTIONS}  tree_defs.o init_tree.o input_output.o out_tree_mod.o show_trees.o -o show_trees.exe

init_tree.o : tree_defs.f90

input_output.o : init_tree.f90 tree_defs.f90

out_tree_mod.o : tree_defs.f90

show_trees.o : init_tree.f90 input_output.f90 out_tree_mod.o tree_defs.o

