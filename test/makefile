C = icc -g -traceback -gdwarf-2 
FORTRAN = ifort -g -gdwarf-2

MKLDIR = /nas02/apps/intel-15.0/cc/mkl

CFLAGS = $(PROFILE) -heap-array -i8 -openmp -w -fast -DMKL_LP64 -m64 \
	 -I$(MKLDIR)/include/intel64/lp64 \
	 -I$(MKLDIR)/include

LFLAGS = $(PROFILE) -L$(MKLDIR)/lib/intel64 -lmkl_intel_lp64 \
	 -lmkl_intel_thread -lmkl_blas95_lp64 -lmkl_lapack95_lp64 \
	 -liomp5 -lmkl_core  -lm -lpthread -traceback

LAPFLAGS = -lcilkrts -L../lib -lhydro_fmm 
EXEC = shelltest1

lap: 
	$(C) $(CFLAGS) -c -I../include $(EXEC).c
	$(C) $(EXEC).o $(LAPFLAGS) -o $(EXEC) $(LFLAGS)

shelltest1:
	./(EXEC)
clean:
	rm -f *.o
	rm -f *.txt
	rm -f $(EXEC)
	rm -f *~


