#You may change the compile option
#MAC
CC = gcc-10 -Os -mcmodel=medium -fopenmp -w
#Linux 1
#CC = gcc -Os -mcmodel=medium -fopenmp -w
#Linux 2
#CC = gcc -Os -mcmodel=medium -std=c99 -fopenmp -w
LIBS = -lm
CFLAGS = -L/usr/local/sac/lib -lsac -lsacio -lm

all: MatchLocate2 SelectFinal SHIFT lsac ccsacc clean

BIN = ../bin

MatchLocate2: MatchLocate2.o sacio.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS) $(CFLAGS)

SelectFinal: SelectFinal.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

SHIFT: SHIFT.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

lsac: lsac.o sacio.o
	$(CC) -o $(BIN)/$@ $^ $(LIBS)

ccsacc: ccsacc.o
	$(CC) -o $(BIN)/$@ $^ $(CFLAGS)

clean:
	rm -f *.o
