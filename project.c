#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <stdbool.h> 

struct process{
	int pid; //process ID
	int AT; //arrival time
	int WT; //Waittime
	int BT; //burst time
	int preemptions;
	int TAT;
};


int main(int argc,char** argv) 
{ 

    if (argc < 7){
        fprintf(stderr, "not enough input arguments\n.");
        return EXIT_FAILURE;
    }
    else if (argc > 8){
        fprintf(stderr, "too many input arguments\n.");
        return EXIT_FAILURE;
    }

    /* 
    first input argument is the see for random number generator
        use srand48() before each scheduling algorithm
        use drand48() obtain the next value in range[0.0, 1.0)
    */
    int rand_num = atoi(argv[1]);

    /*
    second input argument is the lambda for calculating an exponential distribution of interarrival times
    third input argument is the upper bound for exponential distribution numbers
        use in the ceiling function 
    */
    double lambda = atof(argv[2]);
    double tail = atof(argv[3]);
    
    /*
    fourth input argument is the number of processes to simulate
    Process IDs are A->Z, so 26 at most
    */
    int num_of_proc = atoi(argv[4]);
    if (num_of_proc > 26){
        fprintf(stderr, "too many processes\n.");
        return EXIT_FAILURE;
    }

    /*
    fifth input argument is t_cs in milliseconds to perform a context switch
        the 1st half of t_cs is the time to remove given process from CPU
        the 2nd half of t_cs is the time to bring the next process to CPU
        positive even integer
    */
    int context_switch = atoi(argv[5]);
    if (context_switch <= 0 || context_switch % 2 == 1){
        fprintf(stderr, "invalid context_switch, need a positive even integer\n.");
        return EXIT_FAILURE;
    }

    /*
    sixth input agurment is the constant alpha for exponential averaging
        use ceiling function when calculating  τ
    */
    double alpha = atof(argv[6]);

    /*
    seventh input argument is the t_slice for time slice in RR in milliseconds
    eighth input argument is rr_add the flag for whether processes are added to the end or beginning of the ready queue
    */
    int time_slice = atoi(argv[7]);
    
    char *begin_or_end = "END";
    begin_or_end = argv[8];

    if ( !(begin_or_end == "BEGINNING" || begin_or_end == "END" || begin_or_end == NULL) ){
        fprintf(stderr, "invalid begin_or_end: %s\n.", begin_or_end);
        return EXIT_FAILURE;
    }

    // ==== finish parcing the input parameters ====

    // start the random number generator

    // the following code is for each process ->
    
    srand48(rand_num);
    bool valid = false;
    float t_arrive; int num_CPU_burst;


    while (valid == false){
        t_arrive = floor( -log( drand48() ) / lambda );
        if (t_arrive < tail){
            valid = true;
        }
    }
    printf("initial process arrival time: %f\n", t_arrive);

    valid = false;
    while (valid == false){
        num_CPU_burst = (int) (-log( drand48() ) / lambda * 100) + 1;
        if (num_CPU_burst < tail){
            valid = true;
        }
    }


    printf("number of CPU bursts: %d\n", num_CPU_burst);

    int t_CPU_burst; int t_IO_burst;
    for ( int i = 0; i < num_CPU_burst; i++){

        if ( i == num_CPU_burst -1 ){
            valid = false;
            while (valid == false){
                t_CPU_burst = ceil( -log( drand48() ) / lambda  );
                if (t_CPU_burst < tail){
                    valid = true;
                }
            }
            printf("%d (last) actual CPU burst %d, actual IO busrt %d\n", i+1, t_CPU_burst, t_IO_burst);
            break;
        }
        
        valid = false;
        while (valid == false){
            t_CPU_burst = ceil( -log( drand48() ) / lambda  );
            if (t_CPU_burst < tail){
                valid = true;
            }
        }
        
        valid = false;
        while (valid == false){
            t_IO_burst = ceil( -log( drand48() ) / lambda  );
            if (t_IO_burst < tail){
                valid = true;
            }
        }


        printf("%d actual CPU burst %d, actual IO busrt %d\n", i+1, t_CPU_burst, t_IO_burst);
        

    }


    return 0;
}
