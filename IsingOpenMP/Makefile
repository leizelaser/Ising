### module load mpi/openmpi/3.1.2-gcc 
### CC = mpiCC
### module load mpi/mpich/3.1.4-gcc8.1 mpicxx
CC  = 	g++
SRC =  test_omp.cc
OFLAGS = -fopenmp -lpthread
EXE = MC_Ising
default:
	$(CC) -O3 $(SRC) -o $(EXE) -lm $(OFLAGS)
debug:
	$(CC) -g $(SRC) -o $(EXE) -lm $(OFLAGS)
singlecore:
	$(CC) -O3 $(SRC) -o $(EXE) -lm
singlecore-debug:
	$(CC) -g $(SRC) -o $(EXE) -lm
