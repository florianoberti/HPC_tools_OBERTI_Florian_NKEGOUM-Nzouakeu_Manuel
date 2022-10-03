CC = gcc
CFLAGS = -lm -lopenblas -llapacke
EXEC = monProg
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all:	$(EXEC)

monProg:	$(OBJ)
		$(CC) -o $@ $^ $(CFLAGS)


%.o: %.c
	$(CC) -o $@ -c  $< $(CFLAGS)
clean: 
	rm -f *.o core

mrproper:	clean
		rm -f $(EXEC)
