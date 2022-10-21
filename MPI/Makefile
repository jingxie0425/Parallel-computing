CC = gcc
MPICC = mpicc
NVCC = nvcc
FLAGS = -std=c99 -Wall -O3
CPP = g++
MPICXX = mpic++
FLAGSCXX = -std=c++14

1d:
	$(MPICXX) $(FLAGSCXX) ass1_mpi.cpp

2d:
	$(MPICC) $(FLAGS) ass2_mpi.c -lm

inputgen:
	$(CPP) $(FLAGSCXX) 2d_input_gen.cpp -lm

clean:
	rm -f *.out
