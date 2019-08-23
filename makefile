all : src/parallel_gmres.f90 test/test_gmres.f90
	mpif90 -O3 -fPIC -c src/parallel_gmres.f90 -o obj/pgmres.o
	mpif90 -O3 -fPIC -c test/test_gmres.f90
	mpif90 obj/pgmres.o test_gmres.o -o test/test.out
	mv pgmres.mod obj
	rm *.o
