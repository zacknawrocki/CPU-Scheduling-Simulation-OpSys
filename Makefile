.PHONY: project clean
SOURCES := $(wildcard *.c)

project: 
	gcc -Wall -o project -ggdb3 -pg -O0 project.c settings.c process.c queue.c round_robin.c -fdiagnostics-color=always -lm 

clean: 
	rm -f *.o project a.out
