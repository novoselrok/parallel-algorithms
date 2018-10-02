all: chapel openmp mpi

chapel:
	chpl --fast -o build/main-chapel chapel/main.chpl

openmp:
	gcc-8 -O3 -fopenmp -o build/main-openmp openmp/main.c

mpi:
	mpicc -O3 -o build/main-mpi mpi/main.c

clean:
	rm build/main-chapel build/main-openmp build/main-mpi

.PHONY: chapel openmp mpi
