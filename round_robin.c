#include "round_robin.h"

void FCFS(settings *config) {
    config->rr_queue_push_end = true;
    config->time_slice = INT_MAX;
    RR(config);
}

void RR_transition_running_process(settings *config, int t) {
    queue *q = config->q;
    process** procs = config->procs;
    process *current_proc = queue_peek(q);
    if (current_proc == NULL) return;
    int num_procs = config->num_procs;
    int half_t_cx = config->t_cx;

    int* current_burst = current_proc->burst[current_proc->counter_cpu_burst];
    burst_type current_burst_type = current_proc->current_burst_type;
    int current_cpu_burst_length = current_burst[0];
    int current_io_burst_length = current_burst[1];
    int time_in_burst = t - current_proc->current_burst_start; 
    int time_left_in_burst = current_cpu_burst_length - time_in_burst;
    int time_until_preempt = config->time_slice - time_in_burst;

    if (current_burst_type == CX_OFF && time_in_burst >= half_t_cx) { // if the context was switching off this proc
        // 1a) see if the running process finishes its cx_off switch (either take them off the queue or requeue)
        queue_pop(q);
        if (current_cpu_burst_length == 0 && current_io_burst_length < 0) {
            // this was the last burst, so there is no IO burst and we're done
            current_proc->current_burst_type = IO_BURST;
            current_proc->state = FINISHED;
        } else if (current_cpu_burst_length > 0) {
            // process had been preempted: add back to queue
            queue_push(q, current_proc, config->rr_queue_push_end); // push current proc to back of queue
            current_proc->current_burst_type = CPU_BURST;
            current_proc->state = READY;
            current_proc->current_burst_start = -1;
        } else { // if current_cpu_burst == 0 and io burst isn't negative
            // process had finished its cpu burst -- now do IO
            current_proc->current_burst_type = IO_BURST;
            current_proc->state = BLOCKED;
            current_proc->current_burst_start = t;
        }

    } else if (current_burst_type == CX_ON && time_in_burst >= half_t_cx) {
        // 1b) see if the process using the cpu finishes its cx_on switch (mark that the cpu burst has begun)
        current_proc->current_burst_type = CPU_BURST;
        current_proc->current_burst_start = t;
        // state is already RUNNING

    } else if (current_burst_type == CPU_BURST && time_left_in_burst <= 0) {
        // 2) see if the process using the cpu finishes its burst (if so, context switch them out)
        current_proc->current_burst_type = CX_OFF;
        current_proc->current_burst_start = t;
        current_burst[0] = 0;

    } else if (queue_length(q) > 1 && time_until_preempt <= 0 && current_burst_type == CPU_BURST) {
        // 3) see if any time slices expire (RR mode only); move them to the end of the queue
        current_burst[0] = time_left_in_burst;
        current_proc->current_burst_type = CX_OFF;
        current_proc->current_burst_start = t;
    }
}

void RR_handle_io_done(settngs *config, process *proc, int t) {
    if (proc->state != BLOCKED || proc->current_burst_type != IO_BURST) continue;
    int time_in_burst = t - current_proc->current_burst_start; 
    int current_burst_length = current_proc->burst[current_proc->counter_cpu_burst][0];
    if (time_in_burst < current_burst_length) continue;

    // If we get here, the process was blocked due to being in an IO_BURST, which is now done
    proc->state = READY;
    current_proc->current_burst_start = -1;
    queue_push(q, proc, rr_queue_push_end);
}

void RR_clock_tick(settings *config, int t) {
    process** procs = config->procs;
    int num_procs = config->num_procs;
    int t_cx = config->t_cx;
    int time_slice = config->time_slice;
    bool rr_queue_push_end = config->rr_queue_push_end;
    queue *q = config->q;

    process *current_proc = queue_peek(q);
    // 1a) see if the running process finishes its cx_off switch (either take them off the queue or requeue)
    // 1b) see if the process using the cpu finishes its cx_on switch (mark that the cpu burst has begun)
    // 2) see if the processes using the cpu finishes its burst (either take them off the queue or requeue)
    // 3) see if any time slices expire (RR mode only); move them to the end of the queue
    if (current_proc != NULL) RR_transition_running_process(config, t);

    // if the current proc was popped off queue because it CX_OFFed, start the context switch on the next one
    if (current_proc != queue_peek(q)) {
        current_proc = queue_peek(q);
        current_proc->state = RUNNING;
        current_proc->current_burst_type = CX_ON;
        current_proc->current_burst_start = t;
    }


    // 4) see if any processes finish their IO (requeue them)
    for (int i = 0; i < num_procs; ++i) {
        process *proc = procs[i];
    }

    // 5) see if any processes arrive (add them to queue after breaking ties)
    for (int i = 0; i < num_procs; ++i) {
        process *proc = procs[i];
        if (proc->state != NOT_YET_ARRIVED) continue;
        if (proc->t_arrive > t) continue;

        // If we get here, process should arrive now
        proc->state = READY;
        proc->counter_cpu_burst = 0;
        queue_push(q, proc, rr_queue_push_end);
    }

    // 6) see if there are no more processes (exit simulation)
    // No code needed. Checked in for loop
}

void RR(settings *config) {
    print_config(config);
    return;
    // To start, each process should be in NOT_YET_ARRIVED state
    int t = 0;
    for (t = 0; num_running_procs(config) > 0; ++t) {
        RR_clock_tick(config, t);
    }
}
