.PHONY: project project-no-algo project-fcfs-only clean

C := gcc
CFLAGS := -Wall -fdiagnostics-color=always 
CFLAGS_DEBUG := -ggdb3 -pg -O0
LDFLAGS := -lm

SOURCES := $(wildcard *.c)

project: 
	$(C) $(CFLAGS) $(CFLAGS_DEBUG) -o project $(SOURCES) $(LDFLAGS)

project-no-algo: 
	$(C) $(CFLAGS) $(CFLAGS_DEBUG) -DNO_ALGOS -o project $(SOURCES) $(LDFLAGS)

project-rr-only: 
	$(C) $(CFLAGS) $(CFLAGS_DEBUG) -DNO_SRT -DNO_SJF -DNO_FCFS -o project $(SOURCES) $(LDFLAGS)

clean: 
	rm -f *.o project a.out
