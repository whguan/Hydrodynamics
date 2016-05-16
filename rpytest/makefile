C = icc -O3 #-g 
FORTRAN = ifort

MKLDIR = /nas02/apps/intel-15.0/cc/mkl

CFLAGS = $(PROFILE) 
LFLAGS = $(PROFILE)


LAPFLAGS = -lcilkrts -L../lib -lhydro_fmm 
EXEC = rpytest

lap: 
	$(C) $(CFLAGS) -c -I../include $(EXEC).c
	$(C) $(LFLAGS) $(EXEC).o $(LAPFLAGS) -o $(EXEC)

rpytest:
	./(EXEC)
clean:
	rm -f *.o
	rm -f *.txt
	rm -f $(EXEC)
	rm -f *~


