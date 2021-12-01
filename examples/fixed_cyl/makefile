EXAMPLE = fixed_cyl
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
		 -g -O0 -Wall -Wconversion -fimplicit-none -fbacktrace
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
	FFLAGS = -I$(MKLROOT)/include/intel64/lp64 \
		 -I"${MKLROOT}/include" \
		 -g -warn all -check all -traceback -CB -debug all
#		 -init=snan -O2

# Dynamic linking
#	LFLAGS = ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
#		-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed \
#	 	-lmkl_lp64 -lmkl_sequential -lmkl_core \
#	 	-lpthread -lm -ldl \
#	 	-g -check all -warn all -traceback
#	 	-O3

# Static linking
	LFLAGS = ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
		 -Wl,--start-group \
		 ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
		 ${MKLROOT}/lib/intel64/libmkl_sequential.a \
		 ${MKLROOT}/lib/intel64/libmkl_core.a \
		 -Wl,--end-group -lpthread -lm -ldl \
		 -g -check all -warn all -traceback
#		 -O3
endif

$(EXE) : $(OBJDIR)/vectors.o $(OBJDIR)/contours.o \
	$(OBJDIR)/vortex.o $(OBJDIR)/$(EXAMPLE).o
	$(LD) -o $(EXE) $(OBJDIR)/vectors.o $(OBJDIR)/contours.o \
		$(OBJDIR)/vortex.o $(OBJDIR)/$(EXAMPLE).o $(LFLAGS)

#$(MODDIR)/lapack95.mod : $(MKLROOT)/include/lapack.f90
#	$(FC) -fsyntax-only $< $(FFLAGS)

$(OBJDIR)/vectors.o $(MODDIR)/vectors.mod : $(SRCDIR)/vectors.f90
	$(FC) -c -o $(OBJDIR)/vectors.o $< $(FFLAGS)

$(OBJDIR)/contours.o $(MODDIR)/contours.mod : \
	$(SRCDIR)/contours.f90 $(MODDIR)/vectors.mod
	$(FC) -c -o $(OBJDIR)/contours.o $< $(FFLAGS)

$(OBJDIR)/vortex.o $(MODDIR)/vortex.mod : $(SRCDIR)/vortex.f90 \
	$(MODDIR)/vectors.mod $(MODDIR)/contours.mod #$(MODDIR)/lapack95.mod
	$(FC) -c -o $(OBJDIR)/vortex.o $< $(FFLAGS)

$(OBJDIR)/$(EXAMPLE).o : $(EXAMPLE).f90 $(MODDIR)/vectors.mod \
	$(MODDIR)/contours.mod $(MODDIR)/vortex.mod
	$(FC) -c -o $(OBJDIR)/$(EXAMPLE).o $< $(FFLAGS)

.PHONY:clean cleanall
clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod
cleanall:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(EXE)
