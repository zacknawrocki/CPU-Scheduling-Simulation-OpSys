#include <string.h>
#include <stdio.h>
#include "settings.h"
#include "process.h"

process *copy_process(const process* proc_tmpl) {
    process *proc = malloc(sizeof(process));
    memcpy(proc, proc_tmpl, sizeof(process));
    proc->burst = malloc(sizeof(proc_tmpl->burst));
    memcpy(proc->burst, proc_tmpl->burst, sizeof(proc_tmpl->burst));
    int num_bursts = sizeof(proc_tmpl->burst) / (sizeof(int*));
    for (int i = 0; i < num_bursts; ++i) {
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
    int num_bursts = sizeof(proc->burst) / (sizeof(int*) * 2);
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
    printf("Process %c:\n", proc->id);
    printf("Process %c arrives at: %dms:\n", proc->id, proc->t_arrive);
    printf("Process %c CPU bursts: %d\n", proc->id, proc->num_cpu_burst);
    printf("Process %c current Burst: #%d %s (started at: %dms)\n", proc->id, proc->counter_cpu_burst, burst_types[proc->current_burst_type], proc->current_burst_start);
    // TODO(?): List bursts?
    printf("Process %c tau: %d (next: %d)\n", proc->id, proc->tau, proc->next_tau);
    printf("Process %c turn-around time: %d\n", proc->id, proc->TAT);
    printf("Process %c wait time: %d\n", proc->id, proc->WT);
    printf("Process %c current state: %s\n", proc->id, states[proc->state]);
}
