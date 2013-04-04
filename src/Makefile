CC = g++
CFLAGS = -c
LDFLAGS = -L/Users/ktussey/usr/local/lib -lopenbabel -lboost_program_options -I/Users/ktussey/usr/local/include
GFLAGS= -fopenmp

all: spe

spe:  spe.o read_write_files.o spe_functions.o tanimoto.o
	$(CC) $(LDFLAGS) $(GFLAGS) spe.o spe_functions.o tanimoto.o read_write_files.o -o spe

spe.o: spe.cpp
	$(CC) $(LDFLAGS) $(CFLAGS) $(GFLAGS) spe.cpp -o spe.o

spe_functions.o: spe_functions.cpp
	$(CC) $(LDFLAGS) $(CFLAGS) $(GFLAGS) spe_functions.cpp -o spe_functions.o

tanimoto.o: tanimoto.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) $(GFLAGS) tanimoto.cpp -o tanimoto.o
	
read_write_files.o: read_write_files.cpp
	$(CC) $(CFLAGS) $(LDFLAGS) $(GFLAGS) read_write_files.cpp -o read_write_files.o

clean:
	rm -rf spe.o read_write_files.o tanimoto.o spe_functions.o spe

