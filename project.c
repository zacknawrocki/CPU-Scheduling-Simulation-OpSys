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
    int counter_cpu_burst;
	int **burst; 
    int tau;
    int next_tau;  //SRT need it to calculate tau
	int preemptions;
	int TAT;
	int WT;
} pcs;


void SRT(struct process *ptr_pcs, int num_of_proc, int context_switch, double alpha);
void FCFS(struct process *ptr_pcs, int num_of_proc, int context_switch, double alpha);
void output_file(int algorithm, int avg_BT, int avg_WT, int avg_TAT, int context_switches, int preemptions);

int main(int argc,char** argv) 
{ 

    if (argc < 7){
        fprintf(stderr, "not enough input arguments\n.");
        return EXIT_FAILURE;
    }
    else if (argc > 9){
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
    
    // ===========NEED FIX================
    // this code causes segmentation fault
    // ===================================
    /*
    char *begin_or_end = "END";
    begin_or_end = argv[8];
    
    if (strcmp(begin_or_end, "BEGINNING") == 0 || strcmp(begin_or_end, "END") == 0 || strcmp(begin_or_end, NULL) == 0) {
        fprintf(stderr, "invalid begin_or_end: %s\n.", begin_or_end);
        return EXIT_FAILURE;
    }
    */


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
        ptr_pcs -> counter_cpu_burst = num_CPU_burst;

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
    
    //reset the pointer point to the start of the all_process array
    ptr_pcs = all_processes;
    FCFS(ptr_pcs, num_of_proc, context_switch, alpha);
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

void FCFS(struct process *ptr_pcs, int num_of_proc, int context_switch, double alpha) {

    // Dynamically allocate again, so the original data is not modified
    struct process fcfs_all_processes[num_of_proc];
    struct process *fcfs_ptr_pcs = fcfs_all_processes;

    for (int i = 0; i < num_of_proc; i++) {
    	fcfs_ptr_pcs->id = ptr_pcs->id;
        fcfs_ptr_pcs->t_arrive = ptr_pcs->t_arrive;
        fcfs_ptr_pcs->num_cpu_burst = ptr_pcs->num_cpu_burst;
        printf("Process %c [NEW] (arrival time %d ms) %d CPU bursts\n", fcfs_ptr_pcs->id, fcfs_ptr_pcs->t_arrive, fcfs_ptr_pcs->num_cpu_burst);
        int **fcfs_burst;
        fcfs_burst = calloc(fcfs_ptr_pcs -> num_cpu_burst, sizeof(int *));
        for (int j = 0; j < fcfs_ptr_pcs -> num_cpu_burst; j++) {
            fcfs_burst[j] = calloc(2, sizeof(int));
            fcfs_burst[j][0] = ptr_pcs -> burst[j][0];
            fcfs_burst[j][1] = ptr_pcs -> burst[j][1];
        }
        fcfs_ptr_pcs -> burst = fcfs_burst;
        fcfs_ptr_pcs++;
        ptr_pcs++;
    }
    fcfs_ptr_pcs = fcfs_all_processes;
    struct process *queue[num_of_proc];
    for (int i = 0; i < num_of_proc; i++) {
    	queue[i] = fcfs_ptr_pcs;
    	fcfs_ptr_pcs++;
    }



    // Setup for FCFS Simulation
    int time = 0;
    int wait_times[num_of_proc];
    int turn_around_times[num_of_proc];
    int bursts[num_of_proc];
    int burst_time = 0;
    int curr_index = 0;
    int newest_index = 0;
    int context_switches = 0;
    int avg_BT = 0;
    int avg_WT = 0;
    int avg_TAT = 0;
    printf("time 0ms: Simulator started for FCFS [Q <empty>]");


    // FCFS Simulation
    while(1) {
    	// Check if something can arrive to queue at this given time
        fcfs_ptr_pcs = fcfs_all_processes;
    	for (int i = 0; i < num_of_proc; i++) {
    		if (fcfs_ptr_pcs->t_arrive == time) {
    			// Add to queue
    			fcfs_ptr_pcs->WT = time;
    			queue[newest_index] = fcfs_ptr_pcs;
                bursts[newest_index] = fcfs_ptr_pcs->num_cpu_burst; 
    			newest_index += 1;

    			// Print that you have added to queue
    			printf("time %dms: Process %c arrived; added to ready queue [Q ", time, fcfs_ptr_pcs->id);
    			  // If queue is now empty, after starting new process
    			if (curr_index == newest_index) {
    				printf("<empty>]");
    			} else {
    				for (int j = curr_index; j < newest_index; j++) {
    					if (j != newest_index - 1) {
    						printf("%c ", queue[j]->id);
    					} else {
    						printf("%c", queue[j]->id);
    					}
    				}
    				printf("]\n");
    			} 
    		}
            fcfs_ptr_pcs++;
        }


        // If the current process is done, a new process is started and the queue is updated
        if (burst_time <= 0) {
        	// Every process has been taken care of -- done with FCFS simulation
        	if (curr_index == num_of_proc - 1) {
        		printf("time %dms: Simulator ended for FCFS [Q <empty>]", time);
        		break;
        	}
        	// Iterate to next time if nothing is in the queue
        	if (curr_index == newest_index) {
        		continue;
        	}

        	// Simulate "popping a queue" by iterating to next index
        	curr_index += 1;
        	burst_time = queue[curr_index]->num_cpu_burst;
        	// Store the current process' wait time, to find the average later
        	wait_times[curr_index] = time - queue[curr_index]->WT;
        	// Store the previous process' turn around time, to find the average later
        	turn_around_times[curr_index - 1] = time - queue[curr_index - 1]->WT;

        	// Print the old process has terminated
        	printf("time %dms: Process %c terminated [Q ", time, queue[curr_index - 1]->id);
        	// If queue is now empty, after starting new process
        	if (curr_index - 1 == newest_index) {
        		printf("<empty>]");
        	} else {
    			for (int j = curr_index; j < newest_index; j++) {
    				if (j != newest_index - 1) {
    					printf("%c ", queue[j]->id);
    				} else {
    					printf("%c", queue[j]->id);
    				}
    			}
    			printf("]\n");
    		} 

        	// Print that new process is using CPU
        	printf("time %dms: Process %c started using the CPU for %dms burst [Q ", time, queue[curr_index]->id, queue[curr_index]->num_cpu_burst);
        	// If queue is now empty, after starting new process
        	if (curr_index == newest_index) {
        		printf("<empty>]");
        	} else {
    			for (int j = curr_index; j < newest_index; j++) {
    				if (j != newest_index - 1) {
    					printf("%c ", queue[j]->id);
    				} else {
    					printf("%c", queue[j]->id);
    				}
    			}
    			printf("]\n");
    		} 
        }
    	time++;
    	burst_time--;
	}

    // Calculate averages and output in simout.txt
    for (int i = 0; i < num_of_proc; i++) {
        avg_BT += bursts[i];
        avg_WT += wait_times[i];
        avg_TAT += turn_around_times[i]; 
    } 
    output_file(2, avg_BT / num_of_proc, avg_WT / num_of_proc, avg_TAT / num_of_proc, context_switches, 0);

	// Deallocate FCFS memory
    fcfs_ptr_pcs = fcfs_all_processes;
	for (int i = 0; i < num_of_proc; i++) {
        for (int j=0; j < fcfs_ptr_pcs->num_cpu_burst; j++){
            free(fcfs_ptr_pcs->burst[j]);
        }
        free(fcfs_ptr_pcs->burst);
        fcfs_ptr_pcs++;
    }
}



void SRT(struct process *ptr_pcs, int num_of_proc, int context_switch, double alpha){
    //need to dynamically again for all_process within this function so that we won't modify the original data.

    struct process srt_all_processes[num_of_proc];
    struct process *srt_ptr_pcs = srt_all_processes;

    for (int i = 0; i < num_of_proc; i++){
        srt_ptr_pcs -> id = ptr_pcs -> id;
        srt_ptr_pcs -> t_arrive = ptr_pcs -> t_arrive;
        srt_ptr_pcs -> num_cpu_burst = ptr_pcs -> num_cpu_burst;
        srt_ptr_pcs -> counter_cpu_burst = ptr_pcs -> counter_cpu_burst;
        
        printf("Process %c [NEW] (arrival time %d ms) %d CPU bursts\n", srt_ptr_pcs->id, srt_ptr_pcs->t_arrive, srt_ptr_pcs->counter_cpu_burst);

        int **srt_burst;
        srt_burst = calloc(srt_ptr_pcs -> num_cpu_burst, sizeof(int *));

        for ( int j = 0; j < srt_ptr_pcs -> num_cpu_burst; j++){

            srt_burst[j] = calloc(2, sizeof(int));
            srt_burst[j][0] = ptr_pcs -> burst[j][0];
            srt_burst[j][1] = ptr_pcs -> burst[j][1];

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
    
    // =================
    // simulator starts
    // =================

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

        // at the start of each ms, srt_ptr_pcs points to the first process
        srt_ptr_pcs = srt_all_processes;

        // context switch counter
        if (t_cs > 0){
            t_cs --;
        }
        
        // time test print
        //
        //printf("[TIME] t_run: %d, t_cs: %d, cpu burst: %c, io burst: %c, new burst: %d\n", t_run, t_cs, srt_id_pcs_running_cpu, srt_id_pcs_running_io, new_burst);
        //    

        for (int i = 0; i < num_of_proc; i++){

            // ===========================================================
            // ending conditions for single process and multiple processes
            // ===========================================================

            //chech if a process is finished, remove it from the srt_all_processes array
            
            finish_process = false;    
            srt_ptr_tmp = srt_ptr_pcs -> burst;
            for (int j = 0; j < srt_ptr_pcs -> num_cpu_burst; j++){
                if (srt_ptr_tmp[j][0]== 0 && srt_ptr_tmp[j][1]== -1){
                    finish_process = true;
                    break;
                }
            }
            
            if (finish_process == true){
                // this process is finished
                num_finish_process++;
                //if there's only one process, just break simulator and end;
                if (num_of_proc == 1){
                    finish = true;
                    break;
                }
                // if there are more process, remove this process from the array 
                else {

                    // might have memory issues when removing element from dynamically allocated array 

                    for (int j = i; j < num_of_proc - 1; j++){
                        srt_all_processes[j] = srt_all_processes[j+1];
                    }
                }
            }
            //if all processes are finished, break simulator and end;
            if (num_of_proc == num_finish_process){
                finish = true;
                break;
            }

        
            /*
              If a process arrive at current time
              --  is there any running process?
                    yes -- is the running process tau smaller?
                             yes -- continue running, add new process to queue, NO CS
                             no  -- change the running process to new one (add to first of queue), add running process to second of queue, CS
                    no -- change the running process to the new one(add to the first of queue), CS         
            */
            if (srt_ptr_pcs->t_arrive == t_run){
                
                if (srt_id_pcs_running_cpu != '-'){
                    //new process preempt the running process
                    if (srt_ptr_pcs_running_cpu->tau > srt_ptr_pcs->tau){
                        // move each process in the queue to next pos
                        for (int j = num_of_proc - 1; j > 0; j--){
                            srt_ptr_pcs_cpu_queue[j] = srt_ptr_pcs_cpu_queue[j-1];
                        }
                        // add running process to the queue
                        srt_ptr_pcs_cpu_queue[0] = srt_ptr_pcs_running_cpu;
                        srt_num_pcs_cpu_queue++;
                        srt_id_pcs_running_cpu = '-';
                        // move back again 
                        for (int j = num_of_proc - 1; j > 0; j--){
                            srt_ptr_pcs_cpu_queue[j] = srt_ptr_pcs_cpu_queue[j-1];
                        }
                        struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                        srt_ptr_pcs_cpu_queue[0] = srt_ptr_pcs_current;
                        srt_num_pcs_cpu_queue++;
                        srt_id_pcs_running_cpu = srt_ptr_pcs -> id;
                        t_cs = context_switch / 2;
                        new_burst = true;
                    }
                    //continue running process, add new process to the queue (position depends on tau)
                    else {
                        struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                        for (int j = num_of_proc - 1; j > 0; j--){
                            if (srt_ptr_pcs_cpu_queue[j-1]->tau > srt_ptr_pcs->tau){
                                srt_ptr_pcs_cpu_queue[j] = srt_ptr_pcs_cpu_queue[j-1];
                                continue;
                            }
                            if (srt_ptr_pcs_cpu_queue[j-1]->tau < srt_ptr_pcs->tau){
                                srt_ptr_pcs_cpu_queue[j] = srt_ptr_pcs_current;
                                srt_num_pcs_cpu_queue++;
                                //no cs, no new_burst, no change to running cpu
                            }
                        }
                    }

                }
                // if no process is running
                else{
                    srt_id_pcs_running_cpu = srt_ptr_pcs->id;
                    struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                    srt_ptr_pcs_cpu_queue[srt_num_pcs_cpu_queue] = srt_ptr_pcs_current;
                    srt_num_pcs_cpu_queue++;
                    t_cs = context_switch / 2;
                    new_burst = true;
                }


                /*
                //once a process arrive, add to the queue first
                struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                srt_ptr_pcs_cpu_queue[srt_num_pcs_cpu_queue] = srt_ptr_pcs_current;
                srt_num_pcs_cpu_queue++;
                
                srt_ptr_pcs
                t_cs = context_switch / 2;
                //a new burst is coming after the context_switch
                new_burst = true;
                */
  
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
                
                
            }

            // if in context switch, new process arrive and preempted to wait for context switch, do not do anything just continue;
            if (t_cs != 0){
                continue;
            }


            /*
            if we have previously finished a cpu burst,
                we are either doing nothing (running cpu id is '-')
                              doing context switch (handled by line 378-380 & 532-534 )
                              finished context switch (it has to match the running process id)

            */
            if (new_burst == true){
                
                //
                // printf("[test]process %c finished context switch, it's time to set it run\n",srt_id_pcs_running_cpu);
                //

                //if this current process is the running process.
                if (srt_id_pcs_running_cpu == srt_ptr_pcs->id){

                    // once a process arrive and finish context switch
                    // add it to running process and remove from queue.
                    for (int j = 0; j < num_of_proc; j++){
                        if (srt_ptr_pcs_cpu_queue[j] -> id == srt_id_pcs_running_cpu){
                            //remove this process in the queue
                            for (int k = j; k < num_of_proc-1; k++){
                                srt_ptr_pcs_cpu_queue[k] = srt_ptr_pcs_cpu_queue[k+1];
                            }
                            srt_num_pcs_cpu_queue--;
                        }
                    }

                    //seems to be duplicate code, leave it here for now
                    struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                    srt_ptr_pcs_running_cpu = srt_ptr_pcs_current;

                    //check if a process ends, if end, break

                    //
                    // printf("[test] check if process %c has finished -- ",srt_id_pcs_running_cpu);
                    //
                    
                    //find the first available value in srt_ptr_pcs_running_cpu -> burst
                    // if no available value, this process is finishedcheck if p
                    srt_ptr_tmp = srt_ptr_pcs_running_cpu->burst;
                    int tmp = 0;
                    bool finish_pcs = true;
                    for (int m = 0; m < srt_ptr_pcs_running_cpu->num_cpu_burst; m++){
                        
                        // tmp is the next cpu burst
                        if (srt_ptr_tmp[m][0] != 0){
                            tmp = srt_ptr_tmp[m][0];
                            finish_pcs = false;
                            break;
                        }
                        if (srt_ptr_tmp[m][1] == -1){
                            finish_pcs = true;
                            break;
                        }
                    }

                    printf("%d\n", finish_pcs);

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
 
                        // this process is finished
                        num_finish_process++;
                        //if there's only one process, just break simulator and end;
                        if (num_of_proc == num_finish_process){
                            finish = true;
                            break;
                        }
                        // if there are more process, remove this process from the array 
                        else {
                            // might have memory issues when removing element from dynamically allocated array 

                            for (int j = i; j < num_of_proc - 1; j++){
                                srt_all_processes[j] = srt_all_processes[j+1];
                            }
                        }
                        

                        break;
                    }


                    // if there's is more burst to go for this process
                    // tmp is calculated next cpu burst time
                    srt_ptr_pcs_running_cpu -> next_tau = tmp;                    
                    finish_cpu_burst = false;
                    new_burst = false;


                    // ===============================
                    // formatted print the cpu queue:
                    // ===============================
                    if (srt_num_pcs_cpu_queue == 0){
                        char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                        printf("time %dms: Process %c (tau %dms) started using the CPU with %dms burst remaining [Q %s]\n", t_run, srt_id_pcs_running_cpu, srt_ptr_pcs_running_cpu->tau, tmp, cpu_queue);
                    
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
                        printf("time %dms: Process %c (tau %dms) started using the CPU with %dms burst remaining [Q %s]\n", t_run, srt_id_pcs_running_cpu, srt_ptr_pcs_running_cpu->tau, tmp, cpu_queue);
                    
                    }
                    // ===============================
                    // end of formatted print queue
                    // ===============================

                }
            
            }


            /*
                if code goes here
                    there's no new burst happening
                    there's no new process coming
                    we are doing some normal cpu burst (decreasing the remaining burst time)
                    we are doing nothing (running cpu id = '-')
                if a process finish cpu burst here
                    if io is occupied or io queue is not empty
                        add to the end of io queue
                        set new_burst to true(ready to find the next available burst)
            */

            if (new_burst == false){

                //someone is using cpu and not being preempted
                //continue to subtract the cpu
                if (srt_id_pcs_running_cpu != '-'){
                    
                    // current process is the running process
                    if (srt_id_pcs_running_cpu == srt_ptr_pcs->id){

                        int ** srt_ptr_burst = srt_ptr_pcs_running_cpu -> burst;

                        for (int j = 0; j < srt_ptr_pcs_running_cpu->num_cpu_burst; j++){
                            if(srt_ptr_burst[j][0] == 0 && srt_ptr_burst[j][1] == 0){
                                continue;
                            }

                            //this location is the current burst
                            if (srt_ptr_burst[j][0] != 0){
                                srt_ptr_burst[j][0]--;

                                //
                                // printf("[test] process %c current cpu burst remaining: %d\n", srt_id_pcs_running_cpu, srt_ptr_burst[j][0]);
                                //

                                //if cpu burst not finish, break the loop;
                                if (srt_ptr_burst[j][0] != 0){
                                    break;
                                    
                                }

                                //if finish cpu burst, move to io and perform context switch
                                // need to check if io is occupied & if io queue is empty or not

                                if (srt_ptr_burst[j][0] == 0 && srt_ptr_burst[j][1] != 0){

                                    //
                                    // printf("[test] process %c %dth burst have io remaining, cpu: %d, io: %d\n",srt_id_pcs_running_cpu, j, srt_ptr_burst[j][0], srt_ptr_burst[j][1]);
                                    //

                                    srt_ptr_pcs_running_cpu->counter_cpu_burst--;
                                    
                                    if (srt_id_pcs_running_io == '-'){
                                        srt_ptr_pcs_running_io = srt_ptr_pcs_running_cpu;
                                        srt_id_pcs_running_io = srt_id_pcs_running_cpu;
                                    }
                                    else {
                                        //add to the end of io queue

                                        //TODO: need to modify the counter for occupying IO burst
                                        for (int l = 0; l < num_of_proc; l++){
                                            if (srt_ptr_pcs_io_queue[l]==NULL){
                                                srt_ptr_pcs_io_queue[l] = srt_ptr_pcs_running_cpu;
                                                break;
                                            }
                                        }
                                    }

  
                                    // ===============================
                                    // formatted print the cpu queue:
                                    // ===============================
                                    if (srt_num_pcs_cpu_queue == 0){
                                        char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                                        
                                        printf("time %dms: Process %c (tau %dms) completed a CPU burst; %d bursts to go [Q %s]\n",t_run,srt_ptr_pcs_running_cpu->id, srt_ptr_pcs_running_cpu->tau, srt_ptr_pcs_running_cpu->counter_cpu_burst, cpu_queue);
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
                                        
                                        printf("time %dms: Process %c (tau %dms) completed a CPU burst; %d bursts to go [Q %s]\n",t_run,srt_ptr_pcs_running_cpu->id, srt_ptr_pcs_running_cpu->tau, srt_ptr_pcs_running_cpu->counter_cpu_burst, cpu_queue);
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

                        }
                        
                    }

                }

                
                //do io burst here
                if (srt_id_pcs_running_io != '-'){

                    //only do io if this process is the running io process
                    if (srt_id_pcs_running_io == srt_ptr_pcs->id){

                        int ** srt_ptr_burst = srt_ptr_pcs_running_io -> burst;
                        for (int j = 0; j < srt_ptr_pcs_running_io->num_cpu_burst; j++){

                            if (srt_ptr_burst[j][0] == 0 && srt_ptr_burst[j][1] == -1){
                                finish_process = true;
                                break;
                            }

                            //do io burst for the first available io
                            if (srt_ptr_burst[j][1] != 0){
                                srt_ptr_burst[j][1]--;

                                //
                                // printf("[test] process %c current io burst remaining: %d\n", srt_id_pcs_running_io, srt_ptr_burst[j][1]);
                                //

                                if (srt_ptr_burst[j][1] != 0){
                                    break;
                                    
                                }

                                // a process finish both cpu burst and io burst
                                // need to be added back to the cpu queue
                                // 
                                
                                if (srt_ptr_burst[j][1] == 0){

                                    //if there's cpu running, compare with the runnning cpu check if preemption exists
                                    // else, add to the cpu queue

                                    //
                                    // printf("[test] %c finish one cpu burst - running cpu is %c\n", srt_id_pcs_running_io, srt_id_pcs_running_cpu);
                                    //

                                    if (srt_id_pcs_running_cpu != '-'){
                                        //new process preempt the running process
                                        if (srt_ptr_pcs_running_cpu->tau > srt_ptr_pcs->tau){

                                            //
                                            // printf("[test] %c is running process and tau is bigger than process %c, so being preempted\n", srt_id_pcs_running_cpu, srt_ptr_pcs->id);
                                            //

                                            // move each process in the queue to next pos
                                            for (int k = num_of_proc - 1; k > 0; k--){
                                                srt_ptr_pcs_cpu_queue[k] = srt_ptr_pcs_cpu_queue[k-1];
                                            }
                                            // add running process to the queue
                                            srt_ptr_pcs_cpu_queue[0] = srt_ptr_pcs_running_cpu;
                                            srt_num_pcs_cpu_queue++;
                                            srt_id_pcs_running_cpu = '-';
                                            // move back again 
                                            for (int k = num_of_proc - 1; k > 0; k--){
                                                srt_ptr_pcs_cpu_queue[k] = srt_ptr_pcs_cpu_queue[k-1];
                                            }
                                            struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                                            srt_ptr_pcs_cpu_queue[0] = srt_ptr_pcs_current;
                                            srt_num_pcs_cpu_queue++;
                                            srt_id_pcs_running_cpu = srt_ptr_pcs -> id;
                                            t_cs = context_switch / 2;
                                            new_burst = true;
                                        }
                                        //continue running process, add new process to the queue (position depends on tau)
                                        else {
                                            
                                            //
                                            // printf("[test] %c is running process and tau is smaller than process %c, so process %c add to the queue\n", srt_id_pcs_running_cpu, srt_ptr_pcs->id, srt_ptr_pcs->id);
                                            //
                                            
                                            struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                                            for (int j = num_of_proc - 1; j > 0; j--){
                                                if (srt_ptr_pcs_cpu_queue[j-1]->tau > srt_ptr_pcs->tau){
                                                    srt_ptr_pcs_cpu_queue[j] = srt_ptr_pcs_cpu_queue[j-1];
                                                    continue;
                                                }
                                                //find the proper location in queue and assign
                                                if (srt_ptr_pcs_cpu_queue[j-1]->tau < srt_ptr_pcs->tau){
                                                    srt_ptr_pcs_cpu_queue[j] = srt_ptr_pcs_current;
                                                    srt_num_pcs_cpu_queue++;
                                                    break;
                                                    //no cs, no new_burst, no change to running cpu
                                                }
                                            }
                                        }
                                    }
                                    
                                    // if no process is running, this process's next available burst becomes running process, need cs
                                    else{
                                        
                                        srt_id_pcs_running_cpu = srt_ptr_pcs->id;
                                        struct process *srt_ptr_pcs_current = srt_ptr_pcs;
                                        srt_ptr_pcs_cpu_queue[srt_num_pcs_cpu_queue] = srt_ptr_pcs_current;
                                        srt_num_pcs_cpu_queue++;
                                        t_cs = context_switch / 2;
                                        new_burst = true;
                                        
                                        //
                                        // printf("[test] no running process, process %c becomes running process\n", srt_id_pcs_running_cpu);
                                        //
                                    
                                    }

                                    //struct process *srt_ptr_pcs_current = srt_ptr_pcs_running_io;          
                                    //srt_ptr_pcs_cpu_queue[srt_num_pcs_cpu_queue] = srt_ptr_pcs_current;
                                    //srt_num_pcs_cpu_queue++;
                                    
                                    //srt_ptr_pcs_running_cpu->num_cpu_burst--;
                                        

                                    // ===============================
                                    // formatted print the cpu queue:
                                    // ===============================
                                    if (srt_num_pcs_cpu_queue == 0){
                                        char cpu_queue[] = {'<','e','m','p','t','y','>','\0'};
                                        printf("time %dms: Process %c (tau %dms) completed I/O; added to ready queue [Q %s]\n",t_run,srt_id_pcs_running_io, srt_ptr_pcs_running_io->tau, cpu_queue);
                                    
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
                                        printf("time %dms: Process %c (tau %dms) completed I/O; added to ready queue [Q %s]\n",t_run,srt_id_pcs_running_io, srt_ptr_pcs_running_io->tau, cpu_queue);
                                    
                                    }
                                    // ===============================
                                    // end of formatted print queue
                                    // ===============================


                                    srt_ptr_pcs_running_io = NULL;
                                    srt_id_pcs_running_io = '-';
                                    
                                
                                    /*
                                    // use the first one in queue as the next process for CPU burst (done in previous code line 775-818)
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
                                    */

                                

                                    //printf("time %dms: Process %c (tau %dms) started using the CPU with %dms burst remaining [Q <empty>]\n", t_run, srt_ptr_pcs_running_cpu->id,srt_ptr_pcs_running_cpu->tau, tmp);


                                    break;
                                    
                                }



                                
                            }
                            

                        }
                        
                    }

                }

            }

 



            srt_ptr_pcs++;
        
        }

        /*
        if (t_run > 1500){
            break;
        }
        */

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


// Output for text file
// Four possible integers for algorithm paramater
// 0 = SJF
// 1 = SRT
// 2 = FCFS
// 3 = RR
void output_file(int algorithm, int avg_BT, int avg_WT, int avg_TAT, int context_switches, int preemptions) {
    FILE * file;
    file = fopen("simout.txt", "a");
    char* algorithms[4] = {"SJF", "SRT", "FCFS", "RR"}; 
    fprintf(file, "Algorithm %s\n", algorithms[algorithm]);
    fprintf(file, "-- average CPU burst time: %.3d ms\n", avg_BT);
    fprintf(file, "-- average CPU wait time: %.3d ms\n", avg_WT);
    fprintf(file, "-- average CPU turnaround time: %.3d ms\n", avg_TAT);
    fprintf(file, "-- total number of context switches: %d\n", context_switches);
    fprintf(file, "-- total number of preemptions: %d\n", preemptions);
}










