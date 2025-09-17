LIBS2:= Vector2.o functions.o Jacobi.o Pinning.o Random.o Superfluid.o PDE.o
CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: gin

%.o: %.cpp %.h
	$(CC) -c $< ${CFLAGS}

gin: gin.C ${LIBS2}
	$(CC) $^ -o ../$@ ${CFLAGS}

depinning.x: depinning.C ${LIBS2}
	$(CC) $^ -o ../$@ ${CFLAGS}
	
clean:
	rm *.o *~
