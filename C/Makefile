CC = gcc
TARGET = main
OBJS = pfaffian.o macaulay.o read_cmat.o main.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -g -Wall -O3 -o $@ $^ -lflint -lgmp -lm4ri	

%.o: %.c
	$(CC) -c -g -Wall -O3 -o $@ $<

pfaffian.o: pfaffian.c pfaffian.h
macaulay.o: macaulay.c macaulay.h
read_cmat.o: read_cmat.c read_cmat.h
main.o: main.c pfaffian.c pfaffian.h macaulay.c macaulay.h read_cmat.c read_cmat.h

clean:
	rm *.o