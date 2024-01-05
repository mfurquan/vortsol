EXAMPLE = vortex_cyl
SRCDIR = ../../src
MODDIR = ../../include
OBJDIR = ../../build
#MKLROOT = /opt/intel/mkl

EXE = ./run/$(EXAMPLE)

#FC = gfortran
FC = ifort
LD = $(FC)

F95ROOT = /home/mfurquan/.local/opt/mkl/gfortran

ifeq ($(FC),gfortran)
	FFLAGS = -I$(F95ROOT)/include/intel64/lp64 -m64 \
		 -I"${MKLROOT}/include" -ffree-line-length-0 \
		 -finit-real=nan -fcheck=all -fcheck=bounds \
#		 -g -Og -Wall -Wconversion -fimplicit-none -fbacktrace
#		 -Ofast

# Dynamic linking
#	LFLAGS = ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
#		-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed \
#	 	-lmkl_gf_lp64 -lmkl_sequential -lmkl_core \
#	 	-lpthread -lm -ldl
#	 	-g -fbounds-check -Wall -fbacktrace
#	 	-Ofast

# Static linking
	LFLAGS = ${F95ROOT}/lib/intel64/libmkl_lapack95_lp64.a \
		-Wl,--start-group \
	 	${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \
	 	${MKLROOT}/lib/intel64/libmkl_sequential.a \
	 	${MKLROOT}/lib/intel64/libmkl_core.a \
		-Wl,--end-group -lpthread -lm -ldl
#	 	-g -fbounds-check -Wall -fbacktrace
#	 	-Ofast
endif

ifeq ($(FC),ifort)
#	FFLAGS = -I$(MKLROOT)/include/intel64/lp64 \
#		 -I"${MKLROOT}/include" \
#
	FFLAGS = -init=snan \
		-g -warn all -check all -traceback -CB \
		 -fpe0 -debug all
#		-O2

# Dynamic linking
#	LFLAGS = ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
#		-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed \
#	 	-lmkl_lp64 -lmkl_sequential -lmkl_core \
#	 	-lpthread -lm -ldl \
#	 	-g -check all -warn all -traceback
#	 	-O3

# Static linking
	LFLAGS = ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a \
		 ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
		 -Wl,--start-group \
		 ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
		 ${MKLROOT}/lib/intel64/libmkl_sequential.a \
		 ${MKLROOT}/lib/intel64/libmkl_core.a \
		 -Wl,--end-group -lpthread -lm -ldl \
		 -g -check all -warn all -traceback
#		 -O2
endif

$(EXE) : $(OBJDIR)/global_param.o $(OBJDIR)/utility.o \
	$(OBJDIR)/flow_elements.o $(OBJDIR)/solid_boundary.o \
	$(OBJDIR)/flow.o $(OBJDIR)/$(EXAMPLE).o
	$(LD) -o $(EXE) $(OBJDIR)/global_param.o \
		$(OBJDIR)/utility.o $(OBJDIR)/flow_elements.o \
		$(OBJDIR)/solid_boundary.o $(OBJDIR)/flow.o \
		$(OBJDIR)/$(EXAMPLE).o $(LFLAGS)

$(MODDIR)/lapack95.mod : $(MKLROOT)/include/mkl_lapack.f90
	$(FC) -fsyntax-only $< $(FFLAGS)

$(MODDIR)/blas95.mod : $(MKLROOT)/include/mkl_blas.f90
	$(FC) -fsyntax-only $< $(FFLAGS)

$(OBJDIR)/global_param.o $(MODDIR)/global_param.mod : $(SRCDIR)/global_param.f90
	$(FC) -c -o $(OBJDIR)/global_param.o $< $(FFLAGS)

$(OBJDIR)/utility.o $(MODDIR)/utility.mod : $(SRCDIR)/utility.f90 $(MODDIR)/global_param.mod
	$(FC) -c -o $(OBJDIR)/utility.o $< $(FFLAGS)

$(OBJDIR)/flow_elements.o $(MODDIR)/flow_elements.mod : $(SRCDIR)/flow_elements.f90 $(MODDIR)/flow_elements.mod
	$(FC) -c -o $(OBJDIR)/flow_elements.o $< $(FFLAGS)

$(OBJDIR)/solid_boundary.o $(MODDIR)/solid_boundary.mod : \
	$(SRCDIR)/solid_boundary.f90 $(MODDIR)/global_param.mod \
	$(MODDIR)/utility.mod $(MODDIR)/flow_elements.mod
	$(FC) -c -o $(OBJDIR)/solid_boundary.o $< $(FFLAGS)

$(OBJDIR)/flow.o $(MODDIR)/flow.mod : $(SRCDIR)/flow.f90 \
	$(MODDIR)/global_param.mod $(MODDIR)/utility.mod \
	$(MODDIR)/solid_boundary.mod $(MODDIR)/lapack95.mod \
	$(MODDIR)/blas95.mod
	$(FC) -c -o $(OBJDIR)/flow.o $< $(FFLAGS)

$(OBJDIR)/$(EXAMPLE).o : $(EXAMPLE).f90 $(MODDIR)/global_param.mod \
	$(MODDIR)/utility.mod $(MODDIR)/solid_boundary.mod \
	$(MODDIR)/flow.mod
	$(FC) -c -o $(OBJDIR)/$(EXAMPLE).o $< $(FFLAGS)

.PHONY:clean cleanall
clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod
cleanall:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXE)
