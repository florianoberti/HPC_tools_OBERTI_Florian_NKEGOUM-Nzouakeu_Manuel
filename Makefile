EXEC = prog01 prog02 prog03 prog04 prog05 prog06
CC = gcc
LDLIBS =  -lm -lopenblas -llapacke -fopenmp
CFLAGS = -march=native -falign-loops=16
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all :	$(EXEC)
%o:%.c
	$(CC) -o  $@ -c  $(LDLIBS)
	
prog01 :	$(OBJ)
	$(CC)  -O0 $^  -o  $@ $(LDLIBS)
prog02 :	$(OBJ)
	 $(CC) -O1 $(CFLAGS)  $^  -o $@  $(LDLIBS)

prog03 :	$(OBJ)
	 $(CC) -O2 $(CFLAGS)   $^ -o $@  $(LDLIBS)

prog04 :	$(OBJ)
	 $(CC) -O3 -floop-parallelize-all $(CFLAGS)  $^ -o $@  $(LDLIBS)

prog05 :	$(OBJ)
	 $(CC) -Ofast $(CFLAGS)   $^ -o $@  $(LDLIBS) 

prog06 :	$(OBJ)
	 $(CC) -O3  $(CFLAGS)    $^ -o   $@   $(LDLIBS) 



clean: 
	rm -f *.o core

mrproper:	clean
			rm -f $(EXEC)
