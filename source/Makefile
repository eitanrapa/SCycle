# SCycle build.
#
# Requires PETSc (>= 3.14) configured with HDF5 and MUMPS (HYPRE for the AMG/CG solver options).
#   export PETSC_DIR=/path/to/petsc PETSC_ARCH=arch-name   (PETSC_ARCH empty for a --prefix install)
#   make -C source                     build source/main
#   make -C source WERROR=1            treat compiler warnings as errors
#   make -C source DEBUG_MODULES=-DVERBOSE=2   more run-time tracing (VERBOSE=1 prints one line per step)
#   make -C source clean

all: main

DEBUG_MODULES ?= -DVERBOSE=1
WERROR        ?= 0
CFLAGS         = $(DEBUG_MODULES)
CXXFLAGS       = $(CFLAGS) -std=c++14 -Wall -g
ifeq ($(WERROR),1)
CXXFLAGS      += -Werror
endif
# generate header dependencies (*.d) while compiling
CXXPPFLAGS     = -MMD -MP

OBJECTS := domain.o fault.o genFuncs.o \
 odeSolver.o rootFinder.o \
 linearElastic.o \
 dislocationCreep.o diffusionCreep.o dissolutionPrecipitationCreep.o pseudoplasticity.o powerLaw.o \
 heatEquation.o grainSizeEvolution.o \
 spmat.o sbpOps_m_constGrid.o sbpOps_m_varGrid.o \
 odeSolverImex.o odeSolver_WaveEq.o odeSolver_WaveImex.o pressureEq.o interiorFaultLift.o \
 strikeSlip_linearElastic_qd.o strikeSlip_powerLaw_qd.o \
 strikeSlip_linearElastic_fd.o strikeSlip_linearElastic_qd_fd.o strikeSlip_powerLaw_qd_fd.o

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

main: main.o $(OBJECTS)
	${CXXLINKER} $^ -o $@ ${PETSC_LIB}

.PHONY: all clean
clean::
	rm -f *.o *.d main

-include $(OBJECTS:.o=.d) main.d
