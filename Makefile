CC = gcc
MPICC = mpicc
CFLAGS = -O3 -pthread -lm -lpng -fopenmp -std=c99
DEBUG = -DDEBUG
TIME = -DTIME

all: mpi_static mpi_dynamic omp hybrid
all_debug: mpi_static_debug mpi_dynamic_debug omp_debug hybrid_debug

mpi_static: mpi_static.c
	$(MPICC) mpi_static.c $(CFLAGS) -o mpi_static

mpi_dynamic: mpi_dynamic.c
	$(MPICC) mpi_dynamic.c $(CFLAGS) -o mpi_dynamic

omp: omp.c
	$(MPICC) omp.c $(CFLAGS) -o omp

hybrid: hybrid.c 
	$(MPICC) hybrid.c $(CFLAGS) -o hybrid

# debug program
mpi_static_debug: mpi_static.c
	$(MPICC) mpi_static.c $(CFLAGS) -o mpi_static_debug ${DEBUG} ${TIME}

mpi_dynamic_debug: mpi_dynamic.c
	$(MPICC) mpi_dynamic.c $(CFLAGS) -o mpi_dynamic_debug ${DEBUG} ${TIME}

omp_debug: omp.c
	$(MPICC) omp.c $(CFLAGS) -o omp_debug ${DEBUG} ${TIME}

hybrid_debug: hybrid.c 
	$(MPICC) hybrid.c $(CFLAGS) -o hybrid_debug ${DEBUG} ${TIME}   

# record time program
mpi_static_time: mpi_static.c
	$(MPICC) mpi_static.c $(CFLAGS) -o mpi_static_time ${TIME}

mpi_dynamic_time: mpi_dynamic.c
	$(MPICC) mpi_dynamic.c $(CFLAGS) -o mpi_dynamic_time ${TIME}

omp_time: omp.c
	$(MPICC) omp.c $(CFLAGS) -o omp_time ${TIME}

hybrid_time: hybrid.c 
	$(MPICC) hybrid.c $(CFLAGS) -o hybrid_time ${TIME}

clean:
	rm sequential mpi_static hybrid omp mpi_dynamic \
	mpi_static_debug hybrid_debug omp_debug mpi_dynamic_debug\
	mpi_static_time mpi_dynamic_time omp_time hybrid_time
