C = icc
CFLAGS = -O3
LFLAGS = -nofor_main -lcilkrts -L../lib -ladap_laplace
EXEC = example

all: $(EXEC)

$(EXEC): $(EXEC).o 
	$(C) $< -o $(EXEC) $(LFLAGS)

$(EXEC).o: $(EXEC).c
	$(C) $(CFLAGS) -c -I../include -o $@ $<

clean: 
	rm -f *.o
	rm -f $(EXEC)
	rm -f *~








