all : src/parallel_gmres.f90 test/test_gmres.f90
	mkdir -p lib include obj
	mpif90 -O3 -shared -fPIC -c src/parallel_gmres.f90 -o lib/libpgmres.so
	mv pgmres.mod ./include
	mpif90 -O3 -c test/test_gmres.f90 -o obj/test_gmres.o -I./include
	mpif90 -O3 -o test/test.out obj/test_gmres.o -lpgmres
