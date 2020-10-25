all: hw2

hw2: hw2.c
	gcc -Wall -o hw2 hw2.c -lm -O3 -mcmodel=medium

clean: 
	rm -fr *~ hw2