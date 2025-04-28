CC = gcc 
CFLAGS = -O3  
LDFLAGS = -lm #-fopenmp 
DEBUG = -g -fsanitize=address -lefence -Wall -Wextra #$(" ") #

EXECS = q3_part1 q3_part2
all: $(EXECS)

q3_part1: q3_part1.c v_cycle.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(DEBUG)

q3_part2: q3_part2.c v_cycle.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(DEBUG)

v_cycle.o: v_cycle.c
	$(CC) $(CFLAGS) -c $< $(LDFLAGS) $(DEBUG) 

.PHONY: clean
clean:
	rm -f *.o $(EXECS)

