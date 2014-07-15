# Makefile for RAT (Ramses Analysis Tools)
# 
# OPTIONS for AMR2MAP 
# 
# - JEJE : use module jeje.f90 to convert ramses quantities into cgs (-> includes
#	   getting T, density of neutral gas, etc.)
#
# - SMOOTH : uses Yohan's adaptive smoothing to show off
#

F90    = ifort -free -cpp -O2
OPTION = -DCHEMO #-DSMOOTH #-DJEJE


default: all clean

part2map:	part2map.o
	$(F90) part2map.o -o part2map

amr2map:	ramses_info.o parameter_module.o jeje.o amr2map.o
	$(F90) ramses_info.o parameter_module.o jeje.o amr2map.o -o amr2map #_smooth

ccmassgal: ccmassgal.o
	$(F90) ccmassgal.o -o ccmassgal

amr2subvol: amr2subvol.o
	$(F90) amr2subvol.o -o amr2subvol

part2subvol: part2subvol.o
	$(F90) part2subvol.o -o part2subvol

sfr2subvol: sfr2subvol.o
	$(F90) sfr2subvol.o -o sfr2subvol

youngpart2subvol: youngpart2subvol.o
	$(F90) youngpart2subvol.o -o youngpart2subvol

align_gas: ramses_mod.o align_gas.o
	$(F90) ramses_mod.o align_gas.o -o align_gas

align_star: ramses_mod.o align_star.o
	$(F90) ramses_mod.o align_star.o -o align_star

align_dm: ramses_mod.o align_dm.o
	$(F90) ramses_mod.o align_dm.o -o align_dm

align_gas_v0: ramses_mod.o align_gas_v0.o
	$(F90) ramses_mod.o align_gas_v0.o -o align_gas_v0

treebrick2ascii: ramses_mod.o treebrick2ascii.o
	$(F90) ramses_mod.o treebrick2ascii.o -o treebrick2ascii

all: part2map amr2map ccmassgal amr2subvol part2subvol sfr2subvol youngpart2subvol align_gas align_star align_dm align_gas_v0 treebrick2ascii

#############################################################################
%.o:%.f90
	$(F90) $(OPTION) -c $^ -o $@
#############################################################################
clean :
	rm -f *.o *.mod
#############################################################################
