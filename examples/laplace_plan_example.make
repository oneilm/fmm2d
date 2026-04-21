OS = osx

HOST = gcc-openmp
#HOST = gcc
#HOST = intel-openmp
#HOST = intel

PROJECT = int2-laplace-plan-example

# This makefile assumes the static library is already built.
# Run 'make lib' from the repo root before running this makefile.
#
# STATICLIB points to the repo's local static library.
# Adjust if you have installed fmm2d elsewhere.

STATICLIB = ../lib-static/libfmm2d.a

ifeq ($(HOST),gcc)
    FC = gfortran
    FFLAGS = -fPIC -O3 -funroll-loops -march=native -std=legacy -w
    LIBS = -lm
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran
    FFLAGS = -fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy -w
    LIBS = -lm -lgomp
endif

ifeq ($(HOST),intel)
    FC = ifort
    FFLAGS = -O3 -fPIC -march=native
    LIBS = -lm
endif

ifeq ($(HOST),intel-openmp)
    FC = ifort
    FFLAGS = -O3 -fPIC -march=native -qopenmp
    LIBS = -lm
endif

# dlaran is a test utility not included in the main library
EXTRA_OBJS = ../src/common/dlaran.o

.PHONY: all clean

default: all

all: $(EXTRA_OBJS)
	$(FC) $(FFLAGS) -o $(PROJECT) laplace_plan_example.f90 \
	     $(EXTRA_OBJS) $(STATICLIB) $(LIBS)
	./$(PROJECT)

clean:
	rm -f $(PROJECT) fort.13 *.mod
