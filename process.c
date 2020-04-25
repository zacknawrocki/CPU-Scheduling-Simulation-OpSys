#include <string.h>
#include <stdio.h>
#include "settings.h"
#include "process.h"

process *copy_process(const process* proc_tmpl) {
    process *proc = malloc(sizeof(process));
    memcpy(proc, proc_tmpl, sizeof(process));
    proc->num_cpu_burst = proc_tmpl->num_cpu_burst;
    proc->burst = calloc(sizeof(int*), proc->num_cpu_burst);
    memcpy(proc->burst, proc_tmpl->burst, sizeof(int*) * proc->num_cpu_burst);
    for (int i = 0; i < proc->num_cpu_burst; ++i) {
        proc->burst[i] = calloc(sizeof(int), 2);
        memcpy(proc->burst[i], proc_tmpl->burst[i], sizeof(int) * 2);
    }
    proc->counter_cpu_burst = 0;
    proc->current_burst_start = 0;
    proc->state = NOT_YET_ARRIVED;
    proc->TAT = 0;
    proc->WT = 0;
    return proc;
}

void free_process(process* proc) {
    int num_bursts = proc->num_cpu_burst;
    for (int i = 0; i < num_bursts; ++i) {
        free(proc->burst[i]);
    }
    free(proc->burst);
    free(proc);
}

int num_running_procs(const settings *config) {
    int num_procs = config->num_procs;
    int running_procs = 0;
    for (int i = 0; i < num_procs; ++i) {
        process_state state = config->procs[i]->state;
        if (state != NOT_YET_ARRIVED && state != FINISHED) ++running_procs;
    }
    return running_procs;
}

void print_process(const process *proc) {
    const char *burst_types[4] = { "IO_BURST", "CPU_BURST", "CX_ON", "CX_OFF" };
    const char *states[5] = { "NOT_YET_ARRIVED", "READY", "RUNNING", "BLOCKED", "FINISHED" };
    printf("\tProcess %c:\n", proc->id);
    printf("\t\tProcess %c arrives at: %dms\n", proc->id, proc->t_arrive);
    printf("\t\tProcess %c CPU bursts: %d\n", proc->id, proc->num_cpu_burst);
    printf("\t\tProcess %c current Burst: %s (#%d of %d; started at: %dms)\n", proc->id, burst_types[proc->current_burst_type], proc->counter_cpu_burst, proc->num_cpu_burst, proc->current_burst_start);
    printf("\t\tBursts:\n");
    for (int i = 0; i < proc->num_cpu_burst; ++i) {
        printf("\t\t\t%d) CPU: %dms, IO: %dms\n", i, proc->burst[i][0], proc->burst[i][1]); 
    }
    printf("\t\tProcess %c tau: %d (next: %d)\n", proc->id, proc->tau, proc->next_tau);
    printf("\t\tProcess %c turn-around time: %d\n", proc->id, proc->TAT);
    printf("\t\tProcess %c wait time: %d\n", proc->id, proc->WT);
    printf("\t\tProcess %c current state: %s\n", proc->id, states[proc->state]);
}
