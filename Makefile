execs := q2 q3

CC := gcc 
CFLAGS := -O3  
LDFLAGS := -lm -fopenmp 
DEBUG := -g -fsanitize=address -lefence -Wall -Wextra #$(" ") #

all: $(execs)

q2: q2.c conj_grad.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(DEBUG)

q3: q3.c conj_grad.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(DEBUG)

conj_grad.o: conj_grad.c
	$(CC) $(CFLAGS) -c $< $(LDFLAGS) $(DEBUG) 

.PHONY: clean
clean:
	rm -f *.o $(execs)

