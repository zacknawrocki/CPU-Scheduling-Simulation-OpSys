#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <stdbool.h> 
#include <string.h>

struct process{
	char id; 
	int t_arrive; 
	int num_cpu_burst; 
	int **burst; 
	int preemptions;
	int TAT;
} pcs;


void SRT(struct process *ptr_pcs, int num_of_proc, int context_switch, double alpha);

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
    char ID_list[26] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
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

    //use an array of struct to store all the processes' data
    struct process all_processes[num_of_proc];
    struct process *ptr_pcs = NULL;
    ptr_pcs = all_processes;



    for (int i = 0; i < num_of_proc; i++){
        // assign the id for process
        ptr_pcs -> id = ID_list[i];


        while (valid == false){
            t_arrive =  floor( -log( drand48() ) / lambda );
            if (t_arrive < tail){
                valid = true;
            }
        }
        printf("initial process arrival time: %f\n", t_arrive);

        ptr_pcs -> t_arrive = t_arrive;

        valid = false;
        while (valid == false){
            num_CPU_burst = (int) (-log( drand48() ) / lambda * 100) + 1;
            if (num_CPU_burst < tail){
                valid = true;
            }
        }
        printf("number of CPU burmemset(array, -1, sizeof(array[0][0]) * row * count)sts: %d\n", num_CPU_burst);

        ptr_pcs -> num_cpu_burst = num_CPU_burst;

        //burst store all the bursts' information

        /*
        int burst[num_CPU_burst][2];    

        int (*ptr_burst)[num_CPU_burst][2] = &burst;

        memset(burst, 0, sizeof(int)*num_CPU_burst*2);  

    */

        int **burst;
        burst = calloc(num_CPU_burst, sizeof(int *));


        int t_CPU_burst; int t_IO_burst;

        for ( int i = 0; i < num_CPU_burst; i++){

            burst[i] = calloc(2, sizeof(int));

            if ( i == num_CPU_burst -1 ){
                valid = false;
                while (valid == false){
                    t_CPU_burst = ceil( -log( drand48() ) / lambda  );
                    if (t_CPU_burst < tail){
                        valid = true;
                    }
                }
                printf("%d (last) actual CPU burst %d\n",i+1, t_CPU_burst);
                
                burst[i][0] = t_CPU_burst;
                burst[i][1] = -1;

                break;
            }
            else{
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

                burst[i][0] = t_CPU_burst;
                burst[i][1] = t_IO_burst;

            }

        }


        ptr_pcs -> burst = burst;

        //move to the next process
        ptr_pcs++;

    }

    // finish all tptr_bursthe preparations

    /*
        test code to check if the process data are stored correctly
    */

    ptr_pcs = all_processes;

    for (int j=0; j < num_of_proc; j++){
        int num_cpu =  ptr_pcs->num_cpu_burst ;
        printf("each processes:\nid: %c\narrive time: %d\nnumber of bursts: %d\n",ptr_pcs->id, ptr_pcs->t_arrive, ptr_pcs->num_cpu_burst);
    
        int **ptr_b = ptr_pcs->burst;
        

        
        for (int k=0; k < ptr_pcs->num_cpu_burst; k++){
            printf("%d each burst - actual cpu burst: %d\n               actual io burst: %d\n", k, ptr_b[k][0], ptr_b[k][1]);
        }
        
        ptr_pcs++;
    }

    


    
    //============================================    
    //Implement your own function here
    //============================================  
    
    //FCFS();

    ptr_pcs = all_processes;
    SRT(ptr_pcs, num_of_proc, context_switch, alpha);

    //SJF();

    //RR();

    





    //need to free the dynamically allocated memory at the end of the code

    ptr_pcs = all_processes;

    for(int i = 0; i < num_of_proc; i++){

        for (int j=0; j < ptr_pcs->num_cpu_burst; j++){
            free(ptr_pcs->burst[j]);
        }
        free(ptr_pcs->burst);

        ptr_pcs++;

    }



    return 0;
}



void SRT(struct process *ptr_pcs, int num_of_proc, int context_switch, double alpha){
    //need to dynamically again for all_process within this function so that we won't modify the original data.

    struct process srt_all_processes[num_of_proc];
    struct process *srt_ptr_pcs = srt_all_processes;

    
    for (int i = 0; i < num_of_proc; i++){
        srt_ptr_pcs -> id = ptr_pcs -> id;
        srt_ptr_pcs -> t_arrive = ptr_pcs -> t_arrive;
        srt_ptr_pcs -> num_cpu_burst = ptr_pcs -> num_cpu_burst;

        int **srt_burst;
        srt_burst = calloc(srt_ptr_pcs -> num_cpu_burst, sizeof(int *));

        for ( int i = 0; i < srt_ptr_pcs -> num_cpu_burst; i++){

            srt_burst[i] = calloc(2, sizeof(int));

            srt_burst[i][0] = ptr_pcs -> burst[i][0];
            srt_burst[i][1] = ptr_pcs -> burst[i][1];

        }

        srt_ptr_pcs -> burst = srt_burst;

        srt_ptr_pcs++;
        ptr_pcs++;
    }

    // all the data should have been copyed to srt_all_processes















// TODO: Free the dynamically allocated memory
    
}