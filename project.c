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
    int tau;
    int next_tau;  //SRT need it to calculate tau
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
    // =============================================

    srand48(rand_num);
    bool valid = false;
    float t_arrive; int num_CPU_burst;

    //use an array of struct to store all the processes' data
    struct process all_processes[num_of_proc];
    struct process *ptr_pcs = NULL;
    ptr_pcs = all_processes;


    // iterative through each process to generate and store the data
    for (int i = 0; i < num_of_proc; i++){
        // assign the id for process
        ptr_pcs -> id = ID_list[i];


        while (valid == false){
            t_arrive =  floor( -log( drand48() ) / lambda );
            if (t_arrive < tail){
                valid = true;
            }
        }
        // printf("initial process arrival time: %f\n", t_arrive);

        ptr_pcs -> t_arrive = t_arrive;

        valid = false;
        while (valid == false){
            num_CPU_burst = (int) (-log( drand48() ) / lambda * 100) + 1;
            if (num_CPU_burst < tail){
                valid = true;
            }
        }
        // printf("number of CPU burmemset(array, -1, sizeof(array[0][0]) * row * count)sts: %d\n", num_CPU_burst);

        ptr_pcs -> num_cpu_burst = num_CPU_burst;

        //burst store all the bursts' information

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
                // printf("%d (last) actual CPU burst %d\n",i+1, t_CPU_burst);
                
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
                // printf("%d actual CPU burst %d, actual IO busrt %d\n", i+1, t_CPU_burst, t_IO_burst);

                burst[i][0] = t_CPU_burst;
                burst[i][1] = t_IO_burst;

            }

        }


        ptr_pcs -> burst = burst;

        //move to the next process
        ptr_pcs++;

    }

    // finish all preparations



    /*

    //=============================================================
    // test code to check if the process data are stored correctly
    //=============================================================


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

    */


    
    //============================================    
    //Implement your own function here
    //============================================  
    
    //FCFS();

    //reset the pointer point to the start of the all_process array
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
        
        //printf("SRT each processes:\nid: %c\narrive time: %d\nnumber of bursts: %d\n",srt_ptr_pcs->id, srt_ptr_pcs->t_arrive, srt_ptr_pcs->num_cpu_burst);
    
        printf("Process %c [NEW] (arrival time %d ms) %d CPU bursts\n", srt_ptr_pcs->id, srt_ptr_pcs->t_arrive, srt_ptr_pcs->num_cpu_burst);

        int **srt_burst;
        srt_burst = calloc(srt_ptr_pcs -> num_cpu_burst, sizeof(int *));

        for ( int j = 0; j < srt_ptr_pcs -> num_cpu_burst; j++){

            srt_burst[j] = calloc(2, sizeof(int));
            srt_burst[j][0] = ptr_pcs -> burst[j][0];
            srt_burst[j][1] = ptr_pcs -> burst[j][1];

            //printf("SRT %d each burst - actual cpu burst: %d\n               actual io burst: %d\n", j+1, srt_burst[j][0], srt_burst[j][1]);

        }

        srt_ptr_pcs -> burst = srt_burst;

        // initial guess of tau
        srt_ptr_pcs -> tau = 100;

        srt_ptr_pcs++;
        ptr_pcs++;
    }
    

    //============================================================
    // all the data should have been copyed to srt_all_processes
    //============================================================

    //put srt_ptr_pcs at the head of srt_all_processes
    srt_ptr_pcs = srt_all_processes;

    //simulator timing ,start from t_run = 0ms

    int t_run = 0;

    // need two queues -> cpu queue and io queue.
    struct process *srt_ptr_pcs_cpu_queue[num_of_proc];
    struct process *srt_ptr_pcs_io_queue[num_of_proc];

    int **srt_ptr_tmp;

    

    int srt_num_pcs_cpu_queue = 0;
    int srt_num_pcs_io_queue = 0;

    struct process *srt_ptr_pcs_running_cpu;
    struct process *srt_ptr_pcs_running_io;
    
    char srt_id_pcs_running_cpu = '-';
    char srt_id_pcs_running_io = '-';

    int t_cs = 0;

    bool finish = false;
    bool finish_cpu_burst = false;
    bool new_burst = false;
    bool finish_process = false;
    bool finish_srt;

    int num_finish_process = 0;
    int srt_tau;

    
    // TODO: not consider tie for now 
    // TODO: not consider multi thread for now

    while (finish != true){
        if (t_run == 0){
            printf("time 0ms: Simulator started for SRT [Q <empty>]\n");

        }

        /* 
            for each millisecond in the simulator:
                1. check if any of the process arrive
                    if true -> check if any of the process is running
                        if true -> check if need to preempt
                    if false -> run the arrive process
                2. check if any of the process run
                    if true -> substract 1ms from the current running CPU burst
                3. check if any of the process is on IO burst
                    if true -> substract 1ms from the current running IO burst     

        */
        srt_ptr_pcs = srt_all_processes;

        if (t_cs > 0){
            t_cs --;
        }
        
        // printf("t_run: %d, t_cs: %d, cpu burst: %c, io burst: %c, new burst: %d\n", t_run, t_cs, srt_id_pcs_running_cpu, srt_id_pcs_running_io, new_burst);
            

        for (int i = 0; i < num_of_proc; i++){

            //if a process is finished, remove it from the srt_all_processes array
            finish_process = true;    
            srt_ptr_tmp = srt_ptr_pcs -> burst;
            for (int j = 0; j < srt_ptr_pcs -> num_cpu_burst; j++){
                if (srt_ptr_tmp[j][0]!= 0 || srt_ptr_tmp[j][1]!=0){
                    finish_process = false;
                    break;
                }
            }


            if (finish_process == true){
                num_finish_process++;

                if (num_of_proc == 1){
                    finish = true;
                    break;
                }
                else {
                    for (int j = i; j < num_of_proc - 1; j++){
                        srt_all_processes[j] = srt_all_processes[j+1];
                    }
                }
            }

            if (num_of_proc == num_finish_process){
                finish = true;
                break;
            }

            
            
            
            if (srt_ptr_pcs->t_arrive == t_run){
                
                //once a process arrive, add to the queue first
                struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                srt_ptr_pcs_cpu_queue[srt_num_pcs_cpu_queue] = srt_ptr_pcs_current;
                srt_num_pcs_cpu_queue++;
                srt_id_pcs_running_cpu = srt_ptr_pcs->id;

                t_cs = context_switch / 2;
                //a new burst is coming after the context_switch
                new_burst = true;


                
                //printf("time %dms: Process %c (tau %dms) arrived; added to ready queue [Q %c]\n", t_run, srt_ptr_pcs->id, srt_ptr_pcs->tau, srt_ptr_pcs->id);
                
                // ===============================
                // formatted print the cpu queue:
                // ===============================
                if (srt_num_pcs_cpu_queue == 0){
                    char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                    printf("time %dms: Process %c (tau %dms) arrived; added to ready queue [Q %s]\n", t_run, srt_ptr_pcs->id, srt_ptr_pcs->tau, cpu_queue);
                
                }
                else{
                    char cpu_queue[srt_num_pcs_cpu_queue * 2];
                    memset(cpu_queue, '\0', sizeof(cpu_queue));

                    for (int j = 0; j < srt_num_pcs_cpu_queue; j++){
                        cpu_queue[j*2] = srt_ptr_pcs_cpu_queue[j]->id;
                        if (j > 0){
                            cpu_queue[j*2 - 1] = ' ';
                        }
                    }
                    printf("time %dms: Process %c (tau %dms) arrived; added to ready queue [Q %s]\n", t_run, srt_ptr_pcs->id, srt_ptr_pcs->tau, cpu_queue);

                }
                // ===============================
                // end of formatted print queue
                // ===============================
                

                
                break;
                
            }

            if (t_cs != 0){
                continue;
            }

            if (new_burst == true){
                if (srt_id_pcs_running_cpu == srt_ptr_pcs->id){
                    // once a process arrive and finish context switch, add it to running process and remove from queue.
                    
                    for (int j = 0; j < num_of_proc; j++){
                        if (srt_ptr_pcs_cpu_queue[j] -> id == srt_id_pcs_running_cpu){
                            //remove this process in the queue
                            for (int k = j; k < num_of_proc-1; k++){
                                srt_ptr_pcs_cpu_queue[k] = srt_ptr_pcs_cpu_queue[k+1];
                            }
                            srt_num_pcs_cpu_queue--;
                        }
                    }

                    struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                    srt_ptr_pcs_running_cpu = srt_ptr_pcs_current;
                    srt_id_pcs_running_cpu = srt_ptr_pcs -> id;

                    //next_tau is basically t_i for calculating tau
                    srt_ptr_pcs_running_cpu->next_tau = srt_ptr_pcs->burst[0][0];
                    new_burst = false;

                    
                    // ===============================
                    // formatted print the cpu queue:
                    // ===============================
                    if (srt_num_pcs_cpu_queue == 0){
                        char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                        printf("time %dms: Process %c (tau %dms) started using the CPU with %dms burst remaining [Q %s]\n", t_run, srt_ptr_pcs->id,srt_ptr_pcs->tau ,srt_ptr_pcs->burst[0][0], cpu_queue);
                    
                    }
                    else{
                        char cpu_queue[srt_num_pcs_cpu_queue * 2];
                        memset(cpu_queue, '\0', sizeof(cpu_queue));

                        for (int j = 0; j < srt_num_pcs_cpu_queue; j++){
                            cpu_queue[j*2] = srt_ptr_pcs_cpu_queue[j]->id;
                            if (j > 0){
                                cpu_queue[j*2 - 1] = ' ';
                            }
                        }
                        printf("time %dms: Process %c (tau %dms) started using the CPU with %dms burst remaining [Q %s]\n", t_run, srt_ptr_pcs->id,srt_ptr_pcs->tau ,srt_ptr_pcs->burst[0][0], cpu_queue);
                    
                    }
                    // ===============================
                    // end of formatted print queue
                    // ===============================
                    


                }
            
            }


            if (new_burst == false){

                //someone is using cpu
                if (srt_id_pcs_running_cpu != '-'){
                    int ** srt_ptr_burst = srt_ptr_pcs_running_cpu -> burst;

                    for (int j = 0; j < srt_ptr_pcs_running_cpu->num_cpu_burst; j++){
                        if (srt_ptr_burst[j][0] != 0){
                            srt_ptr_burst[j][0]--;
                            
                            //printf("--remaining time for %c on CPU is %d\n", srt_ptr_pcs_running_cpu->id, srt_ptr_burst[j][0]);

                            if (srt_ptr_burst[j][0] != 0){
                                break;
                                
                            }
                        
                        }

                        //if finish cpu burst, move to io and perform context switch
                        if (srt_ptr_burst[j][0] == 0 && srt_ptr_burst[j][1] != 0){

                            srt_ptr_pcs_running_cpu->num_cpu_burst--;
                            srt_ptr_pcs_running_io = srt_ptr_pcs_running_cpu;
                            srt_id_pcs_running_io = srt_id_pcs_running_cpu;


                            
                            // ===============================
                            // formatted print the cpu queue:
                            // ===============================
                            if (srt_num_pcs_cpu_queue == 0){
                                char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                                
                                printf("time %dms: Process %c (tau %dms) completed a CPU burst; %d bursts to go [Q %s]\n",t_run,srt_ptr_pcs_running_cpu->id, srt_ptr_pcs_running_cpu->tau, srt_ptr_pcs_running_cpu->num_cpu_burst, cpu_queue);
                                // TODO: recalculate tau
                                srt_tau = alpha * (srt_ptr_pcs_running_cpu -> next_tau) + (1 - alpha) * srt_ptr_pcs_running_cpu -> tau;
                                srt_ptr_pcs_running_cpu -> tau = srt_tau;
                                printf("time %dms: Recalculated tau = %dms for process %c [Q %s]\n", t_run, srt_tau, srt_ptr_pcs_running_cpu->id, cpu_queue);
                                printf("time %dms: Process %c switching out of CPU; will block on I/O until time %dms [Q %s]\n", t_run, srt_id_pcs_running_io, t_run + t_cs + srt_ptr_burst[j][1], cpu_queue);

                            }
                            else{
                                char cpu_queue[srt_num_pcs_cpu_queue * 2];
                                memset(cpu_queue, '\0', sizeof(cpu_queue));

                                for (int j = 0; j < srt_num_pcs_cpu_queue; j++){
                                    cpu_queue[j*2] = srt_ptr_pcs_cpu_queue[j]->id;
                                    if (j > 0){
                                        cpu_queue[j*2 - 1] = ' ';
                                    }
                                }
                                
                                printf("time %dms: Process %c (tau %dms) completed a CPU burst; %d bursts to go [Q %s]\n",t_run,srt_ptr_pcs_running_cpu->id, srt_ptr_pcs_running_cpu->tau, srt_ptr_pcs_running_cpu->num_cpu_burst, cpu_queue);
                                // TODO: recalculate tau
                                srt_tau = alpha * (srt_ptr_pcs_running_cpu -> next_tau) + (1 - alpha) * srt_ptr_pcs_running_cpu -> tau;
                                srt_ptr_pcs_running_cpu -> tau = srt_tau;
                                printf("time %dms: Recalculated tau = %dms for process %c [Q %s]\n", t_run, srt_tau, srt_ptr_pcs_running_cpu->id, cpu_queue);
                                printf("time %dms: Process %c switching out of CPU; will block on I/O until time %dms [Q %s]\n", t_run, srt_id_pcs_running_io, t_run + t_cs + srt_ptr_burst[j][1], cpu_queue);

                            }
                            // ===============================
                            // end of formatted print queue
                            // ===============================
                            
                            srt_ptr_pcs_running_cpu = NULL;
                            srt_id_pcs_running_cpu = '-';
                            t_cs = context_switch / 2;
                            
                            /*
                            printf("time %dms: Process %c (tau %dms) completed a CPU burst; %d bursts to go [Q <empty>]\n",t_run,srt_ptr_pcs_running_cpu->id, srt_ptr_pcs_running_cpu->tau, srt_ptr_pcs_running_cpu->num_cpu_burst);
                            // TODO: recalculate tau
                            srt_tau = alpha * (srt_ptr_pcs_running_cpu -> next_tau) + (1 - alpha) * srt_ptr_pcs_running_cpu -> tau;
                            srt_ptr_pcs_running_cpu -> tau = srt_tau;
                            printf("time %dms: Recalculated tau = %dms for process %c [Q <empty>]\n", t_run, srt_tau, srt_ptr_pcs_running_cpu->id);
                            */
              

                            //printf("time %dms: Process %c switching out of CPU; will block on I/O until time %dms [Q <empty>]\n", t_run, srt_id_pcs_running_io, t_run + t_cs + srt_ptr_burst[j][1]);

                            break;
                            
                        }

                    }
                    break;

                }

                if (srt_id_pcs_running_io != '-'){
                    int ** srt_ptr_burst = srt_ptr_pcs_running_io -> burst;
                    for (int j = 0; j < srt_ptr_pcs_running_io->num_cpu_burst; j++){
                        if (srt_ptr_burst[j][1] != 0){
                            srt_ptr_burst[j][1]--;
                            //printf("--remaining time for %c on IO is %d\n", srt_ptr_pcs_running_io->id, srt_ptr_burst[j][1]);


                            if (srt_ptr_burst[j][1] != 0){
                                break;
                                break;
                            }
                            else if (srt_ptr_burst[j][1] == 0){
                                finish_cpu_burst = true;
                            }
                        }

                        if (srt_ptr_burst[j][1] == 0 && finish_cpu_burst == true){

                            struct process *srt_ptr_pcs_current = srt_ptr_pcs_running_io;          
                            srt_ptr_pcs_cpu_queue[srt_num_pcs_cpu_queue] = srt_ptr_pcs_current;
                            srt_num_pcs_cpu_queue++;
                            
                            

                            // ===============================
                            // formatted print the cpu queue:
                            // ===============================
                            if (srt_num_pcs_cpu_queue == 0){
                                char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                                printf("time %dms: Process %c (tau %dms) completed I/O; added to ready queue [Q %s]\n",t_run,srt_ptr_pcs_running_io->id,srt_ptr_pcs_running_io->tau, cpu_queue);
                            
                            }
                            else{
                                char cpu_queue[srt_num_pcs_cpu_queue * 2];
                                memset(cpu_queue, '\0', sizeof(cpu_queue));

                                for (int j = 0; j < srt_num_pcs_cpu_queue; j++){
                                    cpu_queue[j*2] = srt_ptr_pcs_cpu_queue[j]->id;
                                    if (j > 0){
                                        cpu_queue[j*2 - 1] = ' ';
                                    }
                                }
                                printf("time %dms: Process %c (tau %dms) completed I/O; added to ready queue [Q %s]\n",t_run,srt_ptr_pcs_running_io->id,srt_ptr_pcs_running_io->tau, cpu_queue);
                            
                            }
                            // ===============================
                            // end of formatted print queue
                            // ===============================

                            srt_ptr_pcs_running_io = NULL;
                            srt_id_pcs_running_io = '-';
                            
                            
                            //printf("time %dms: Process %c (tau %dms) completed I/O; added to ready queue [Q A]\n",t_run,srt_ptr_pcs_running_io->id,srt_ptr_pcs_running_io->tau);
                            



                            // use the first one in queue as the next process for CPU burst
                            srt_ptr_pcs_running_cpu = srt_ptr_pcs_cpu_queue[0];
                            srt_id_pcs_running_cpu = srt_ptr_pcs_cpu_queue[0] -> id;

                            for (int k = 0; k < num_of_proc; k++){
                                if (srt_ptr_pcs_cpu_queue[k] -> id == srt_id_pcs_running_cpu){
                                    //remove this process in the queue
                                    for (int l = k; l < num_of_proc - 1; l++){
                                        srt_ptr_pcs_cpu_queue[l] = srt_ptr_pcs_cpu_queue[l+1];
                                    }
                                    srt_num_pcs_cpu_queue--;
                                    break;
                                }
                            }

                            t_cs = context_switch / 2;

                            //find the first available value in srt_ptr_pcs_running_cpu -> burst
                            srt_ptr_tmp = srt_ptr_pcs_running_cpu->burst;
                            int tmp = 0;
                            bool finish_pcs = true;
                            for (int m = 0; m < srt_ptr_pcs_running_cpu->num_cpu_burst; m++){
                                if (srt_ptr_tmp[m][0] != 0){
                                    tmp = srt_ptr_tmp[m][0];
                                    finish_pcs = false;
                                    break;
                                }
                            }

                            if (finish_pcs == true){

                                // ===============================
                                // formatted print the cpu queue:
                                // ===============================
                                if (srt_num_pcs_cpu_queue == 0){
                                    char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                                    printf("time %dms: Process %c terminated [Q %s]\n", t_run, srt_ptr_pcs_running_cpu->id, cpu_queue);
                                
                                }
                                else{
                                    char cpu_queue[srt_num_pcs_cpu_queue * 2];
                                    memset(cpu_queue, '\0', sizeof(cpu_queue));

                                    for (int j = 0; j < srt_num_pcs_cpu_queue; j++){
                                        cpu_queue[j*2] = srt_ptr_pcs_cpu_queue[j]->id;
                                        if (j > 0){
                                            cpu_queue[j*2 - 1] = ' ';
                                        }
                                    }
                                    printf("time %dms: Process %c terminated [Q %s]\n", t_run, srt_ptr_pcs_running_cpu->id, cpu_queue);
                                
                                }
                                // ===============================
                                // end of formatted print queue
                                // ===============================

                                

                                break;
                            }

                            srt_ptr_pcs_running_cpu -> next_tau = tmp;
                            new_burst = false;
                            finish_cpu_burst = false;


                            // ===============================
                            // formatted print the cpu queue:
                            // ===============================
                            if (srt_num_pcs_cpu_queue == 0){
                                char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                                printf("time %dms: Process %c (tau %dms) started using the CPU with %dms burst remaining [Q %s]\n", t_run, srt_ptr_pcs_running_cpu->id,srt_ptr_pcs_running_cpu->tau, tmp, cpu_queue);
                            
                            }
                            else{
                                char cpu_queue[srt_num_pcs_cpu_queue * 2];
                                memset(cpu_queue, '\0', sizeof(cpu_queue));

                                for (int j = 0; j < srt_num_pcs_cpu_queue; j++){
                                    cpu_queue[j*2] = srt_ptr_pcs_cpu_queue[j]->id;
                                    if (j > 0){
                                        cpu_queue[j*2 - 1] = ' ';
                                    }
                                }
                                printf("time %dms: Process %c (tau %dms) started using the CPU with %dms burst remaining [Q %s]\n", t_run, srt_ptr_pcs_running_cpu->id,srt_ptr_pcs_running_cpu->tau, tmp, cpu_queue);
                            
                            }
                            // ===============================
                            // end of formatted print queue
                            // ===============================


                            //printf("time %dms: Process %c (tau %dms) started using the CPU with %dms burst remaining [Q <empty>]\n", t_run, srt_ptr_pcs_running_cpu->id,srt_ptr_pcs_running_cpu->tau, tmp);


                            break;
                            
                        }

                    }
                    break;

                }

            }

 



            srt_ptr_pcs++;
        
        }




        t_run++;

    }


    printf("time %dms: Simulator ended for SRT [Q <empty>]\n", t_run);

    // Free the dynamically allocated memory

    srt_ptr_pcs = srt_all_processes;

    for(int i = 0; i < num_of_proc; i++){

        for (int j=0; j < srt_ptr_pcs->num_cpu_burst; j++){
            free(srt_ptr_pcs->burst[j]);
        }
        free(srt_ptr_pcs->burst);

        srt_ptr_pcs++;

    }

}