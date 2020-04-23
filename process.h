#pragma once

typedef enum {
	NOT_YET_ARRIVED, READY, RUNNING, BLOCKED, FINISHED
} process_state;

typedef enum {
    IO_BURST, CPU_BURST, CX_ON, CX_OFF
} burst_type;

typedef struct process {
    char id; 
    int t_arrive; 
    int num_cpu_burst; 
    int counter_cpu_burst;
    burst_type current_burst_type;
    int current_burst_start;
    int **burst; 
    int tau;
    int next_tau;  //SRT need it to calculate tau
    //int preemptions;
    int TAT;
    int WT;
    process_state state;
} process;

process *copy_process(const process* proc_tmpl);
void free_process(process* proc);
void print_process(const process *proc);
