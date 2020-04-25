#pragma once
#include <stdio.h>
// Based on https://stackoverflow.com/questions/679979/how-to-make-a-variadic-macro-variable-number-of-arguments
#ifdef DEBUG
    #define DEBUG_PRINT(...) printf(__VA_ARGS__)
#else
    #define DEBUG_PRINT(...) ((void)0)
#endif 

#ifdef NO_ALGOS
    #ifndef NO_RR
        #define NO_RR
    #endif
    #ifndef NO_FCFS
        #define NO_FCFS
    #endif
    #ifndef NO_SJF
        #define NO_SJF
    #endif
    #ifndef NO_SRT
        #define NO_SRT
    #endif
#endif
