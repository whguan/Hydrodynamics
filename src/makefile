CC = icc 
FC = ifort 

CFLAGS = $(PROFILE) -g -I./include  -gdwarf-2#-O3
FFLAGS = $(PROFILE) -g -gdwarf-2# -O3
LFLAGS = $(PROFILE)

COMMON = ./common
OBJDIR = ./build
LAPDIR = ./laplace
RPYDIR = ./rpy
HYDRO  = ./hydro

COMC = $(wildcard $(COMMON)/*.c)
LAPC = $(wildcard $(LAPDIR)/*.c)
RPYC = $(wildcard $(RPYDIR)/*.c)
HYDC = $(wildcard $(HYDRO)/*.c)
HYDF = $(wildcard $(HYDRO)/*.f)

COMOBJ  = $(patsubst $(COMMON)/%.c, $(OBJDIR)/%.o, $(COMC))
LAPCOBJ = $(patsubst $(LAPDIR)/%.c, $(OBJDIR)/%.o, $(LAPC))
RPYOBJ  = $(patsubst $(RPYDIR)/%.c, $(OBJDIR)/%.o, $(RPYC))
HYDROBJ = $(patsubst $(HYDRO)/%.c, $(OBJDIR)/%.o, $(HYDC))
HYDOBJF = $(patsubst $(HYDRO)/%.f, $(OBJDIR)/%.o, $(HYDF))

LIB_LAP = ./build/libhydro_fmm.a

lib_lap: $(COMOBJ) $(LAPCOBJ) $(RPYOBJ) $(HYDROBJ) $(HYDOBJF)
	ar -cvq $(LIB_LAP) $(COMOBJ) $(RPYOBJ) $(LAPCOBJ) $(HYDROBJ) $(HYDOBJF)


$(OBJDIR)/%.o: $(COMMON)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

$(OBJDIR)/%.o: $(LAPDIR)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

$(OBJDIR)/%.o: $(RPYDIR)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

$(OBJDIR)/%.o: $(HYDRO)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<

$(OBJDIR)/%.o: $(HYDRO)/%.f
	$(FC) -c $(FFLAGS) -o $@ $<

clean:
	rm -f $(OBJDIR)/*.o
	#rm -f $(LIB_LAP) $(LIB_YUK)
	rm -f $(LIB_LAP)
#rm -f *~ $(LAPDIR)/*~ $(YUKDIR)/*~ $(COMMON)/*~ $(HYDRO)/*~
	rm -f *~ $(LAPDIR)/*~ $(COMMON)/*~ $(HYDRO)/*~
	rm -f ./include/*~


