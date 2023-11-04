# CC and CFLAGS are varilables
CC = g++
CFLAGS = -c
AR = ar
ARFLAGS = rcv
# -c option ask g++ to compile the source files, but do not link.
# -g option is for debugging version
# -O2 option is for optimized version
DBGFLAGS = -g -D_DEBUG_ON_
OPTFLAGS = -O2

all	: bin/NTU_MPS
	@echo -n ""

# optimized version
bin/NTU_MPS	: main.o
			$(CC) $(OPTFLAGS) main.o -Llib -o bin/mps
main.o 	: src/main.cpp 
			$(CC) $(CFLAGS) $< -Ilib -o $@
MPS_tool.o	: src/MPS_tool.cpp src/MPS_tool.h
			$(CC) $(CFLAGS) $(OPTFLAGS) $< -o $@

# clean all the .o and executable files
clean:
		rm -rf *.o lib/*.a lib/*.o bin/*

