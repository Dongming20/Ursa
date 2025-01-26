F90?=ifx

help:
	@echo "----------------------------------------------------------------------------"
	@echo "Usage: make F90=<f90> {configure,compute,all}"
	@echo ""
	@echo "Example: 'make all' will create all the Ursa executable files in ../bin"
	@echo ""
	@echo "Other possible options:"
	@echo "--- <f90> is your own Fortran90 compiler (possible choices: ifx, gfortran). Default is ifx"
	@echo ""
	@echo "REMARK: -Ursa must be linked with the FEAST eigenvalue library: www.feast-solver.org"
	@echo "          -FEAST must be compiled with the same <f90> and <mpi> options."
	@echo "        -Ursa also needs the Tetgen library (GNU-GPL License): wias-berlin.de/software/tetgen"
	@echo "          -install the executable tetgen in your PATH or copy it in the /bin Ursa directory."
	@echo "----------------------------------------------------------------------------"


NAME_compute = ursa_compute
NAME_configure = ursa_configure


configure: $(NAME_configure)

compute: $(NAME_compute)

all: configure compute clean


OBJS_compute = class_linkedlist.o tools.o basisfunctions.o potentials.o potentials_gw.o ursa_compute.o
OBJS_configure = tools.o ursa_configure.o




## for FEAST
F90FLAGS1 = -L${FEASTROOT}/lib/x64 -lfeast

# intel fortran
ifeq ($(F90),ifx)
## intel MKL
F90FLAGS2 = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl -I"${MKLROOT}/include"
F90FLAGS2 += -O3
## Libxc
#F90FLAGS3 = -I${includedir} -O3 -qopenmp
#F90FLAGS4 = -L${libdir} -lxcf90 -lxc
#F90FLAGS2 += -g -traceback -check all -check bounds -check uninit -ftrapuv -debug all -gen-interfaces  #-heap-arrays 1
endif

# gnu fortran
ifeq ($(F90),gfortran)
## intel MKL
F90FLAGS2 = -m64  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -m64  -I"${MKLROOT}/include" -fallow-argument-mismatch
#F90FLAGS2 += -O3 -fopenmp
## Libxc
#F90FLAGS3 = -I${includedir2} -fallow-argument-mismatch
#F90FLAGS4 = -L${libdir2} -lxcf90 -lxc
endif


# .SUFFIXES:
.SUFFIXES: .f90 .o
.PHONY: clean all


%.o: %.f90
	$(F90) -c $< -o $@


$(NAME_compute) : $(OBJS_compute)
	$(F90) -o $(NAME_compute) $(OBJS_compute) $(F90FLAGS1) $(F90FLAGS2)
	mv $(NAME_compute) ../bin
	
$(NAME_configure) : $(OBJS_configure)
	$(F90) -o $(NAME_configure) $(OBJS_configure) $(F90FLAGS2)
	mv $(NAME_configure) ../bin


clean:
	-@rm  *.o *.mod


