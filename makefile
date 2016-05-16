all: 
	make clean hydro

hydro:
	cd ./src; make lib_lap
	cp ./src/build/libhydro_fmm.a ./lib/
	cd ./test; make clean lap; ./shelltest1  #laptest3 #laptest6
clean:
	rm -f *~
	(cd ./src; make clean)
#(cd ./example; make clean)
	(cd ./test; make clean)
	(cd ./lib; rm -f *.a)

testall:
	make all | grep "\*\*\* "

profilelap:
	make clean
	cd ./src; make lib_lap PROFILE=-pg
	cp ./src/build/libhydro_fmm.a ./lib/
	cd ./test; make clean lap PROFILE=-pg
	rm -f gmon.out
	./test/test_adap_fmm -n 3000000 -d 3 -a 3 -s 30
