
.PHONY: help all build config clean veryclean
.SUFFIXES: .f90 .x

help :
	#
	#    help      : print this help
	#    all       : config then build
	#    config    : configure ./ext directory
	#    build     : build all source files
	#    clean     : remove created files
	#    veryclean : remove ./ext library files
	#

#=================================================================

#LOC_IC	 = /opt/intel
#LOC_IC   = /APP/enhpc/compiler/intel/v11.7/256

FC	= ifort
FCFLAGS = 

#=================================================================

SRC	= vasp_locpot vasp_chg_split vasp_chg_sum

#=================================================================

LOC_EXT	= ./ext
LIB_EXT = -I$(LOC_EXT) -L$(LOC_EXT) -llocfort 

#=================================================================

TARGETS	= $(addprefix ../, $(addsuffix .x, ${SRC}))

#=======================================================================

../%.x : %.f90
	$(FC) $(FCFLAGS) -o $@ $< $(LIB_EXT)

#=================================================================

all : config build

config : 
	make -C $(LOC_EXT) lfort

build : $(TARGETS)


clean :
	rm -f *.o *.mod $(TARGETS)

veryclean: clean
	make -C $(LOC_EXT) clean
