
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

#=================================================================

.PHONY: help all build config clean veryclean
.SUFFIXES: .f90 .x

help :
	#
	#    help      : print this help
	#    all       : make library file '$(LIBNAME).a' for fortran
	#    clean     : remove created files
	#    veryclean : 'clean' with mt_stream created file
	#

all : config build

build : $(TARGETS)

#=======================================================================

config : 
	make -C $(LOC_EXT) lfort

../%.x : %.f90
	$(FC) $(FCFLAGS) -o $@ $< $(LIB_EXT)

clean :
	rm -f *.o *.mod $(TARGETS)

veryclean: clean
	make -C $(LOC_EXT) clean
