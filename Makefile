.PHONY: project debug debug-no-algo clean

C := gcc
CFLAGS := -Wall -fdiagnostics-color=always 
CFLAGS_DEBUG := -ggdb3 -pg -O0 -DDEBUG
LDFLAGS := -lm

SOURCES := $(wildcard *.c)

project: 
	$(C) $(CFLAGS) -o project $(SOURCES) $(LDFLAGS)

debug: 
	$(C) $(CFLAGS) $(CFLAGS_DEBUG) -o project $(SOURCES) $(LDFLAGS)

debug-no-algo: 
	$(C) $(CFLAGS) $(CFLAGS_DEBUG) -DNO_ALGOS -o project $(SOURCES) $(LDFLAGS)

clean: 
	rm -f *.o project a.out
