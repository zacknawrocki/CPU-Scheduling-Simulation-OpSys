#include <stdio.h>
#include "round_robin.h"

void FCFS(settings *config) {
    config->rr_queue_push_end = true;
    config->time_slice = INT_MAX;
    RR(config);
}

void RR_transition_running_process(settings *config, process *proc, int t) {
    queue *q = config->q;
    if (proc == NULL) return;
    int half_t_cx = config->t_cx / 2;

    int* current_burst = proc->burst[proc->counter_cpu_burst];
    burst_type current_burst_type = proc->current_burst_type;
    int current_cpu_burst_length = current_burst[0];
    int current_io_burst_length = current_burst[1];
    int time_in_burst = t - proc->current_burst_start; 
    int time_left_in_burst = current_cpu_burst_length - time_in_burst;
    int time_until_preempt = config->time_slice - time_in_burst;

    if (current_burst_type == CX_OFF && time_in_burst >= half_t_cx) { // if the context was switching off this proc
        // 1a) see if the running process finishes its cx_off switch (either take them off the queue or requeue)
        queue_pop(q);
        if (current_cpu_burst_length > 0) {
            // process is being preempted: add to end of queue
            queue_push(q, proc, config->rr_queue_push_end); // push current proc to back of queue
            proc->current_burst_type = NO_BURST;
            proc->state = READY;
            proc->current_burst_start = -1;
            //printf("* time %dms: Process %c finished context-switching away after being preempted. It has been added back to the end of the queue and marked as READY ", t, proc->id);
            //print_queue_items(config->q);
        } else if (current_io_burst_length < 0) {
            // this was the last burst, so there is no IO burst and we're done
            proc->current_burst_type = NO_BURST;
            proc->state = FINISHED;
            //printf("* time %dms: Process %c finished context-switching away after its final CPU burst. Marking as finished and removed from queue. ", t, proc->id);
            //print_queue_items(config->q);
        } else { // if current_cpu_burst == 0 and io burst isn't negative
            // process had finished its cpu burst -- now do IO
            proc->current_burst_type = IO_BURST;
            proc->state = BLOCKED;
            proc->current_burst_start = t;
            //printf("* time %dms: Process %c finished context-switching away after finishing CPU burst. Will now begin IO burst. ", t, proc->id);
            //print_queue_items(config->q);
        }

    } else if (current_burst_type == CX_ON && time_in_burst >= half_t_cx) {
        // 1b) see if the process using the cpu finishes its cx_on switch (mark that the cpu burst has begun)
        proc->current_burst_type = CPU_BURST;
        proc->current_burst_start = t;
        proc->state = RUNNING;
        if (proc->cpu_burst_preempted) {
            printf("time %dms: Process %c started using the CPU with %dms burst remaining ", t, proc->id, current_cpu_burst_length);
        } else {
            printf("time %dms: Process %c started using the CPU for %dms burst ", t, proc->id, current_cpu_burst_length);
        }
        print_queue_items(config->q);

    } else if (current_burst_type == CPU_BURST && time_left_in_burst <= 0) {
        // 2) see if the process using the cpu finishes its burst (if so, context switch them out)
        proc->current_burst_type = CX_OFF;
        proc->current_burst_start = t;
        proc->cpu_burst_preempted = false;
        current_burst[0] = 0;
        int bursts_left = (proc->num_cpu_burst - proc->counter_cpu_burst) - 1;

        if (current_burst[1] < 0) {
            // this was the last burst, so there is no IO burst and we're done
            proc->state = FINISHED;
            printf("time %dms: Process %c terminated ", t, proc->id);
            print_queue_items(config->q);
        } else {
            printf("time %dms: Process %c completed a CPU burst; %d burst%s to go ", t, proc->id, bursts_left, bursts_left == 1 ? "" : "s");
            print_queue_items(config->q);
            printf("time %dms: Process %c switching out of CPU; will block on I/O until time %dms ", t, proc->id, current_io_burst_length + half_t_cx + t);
            print_queue_items(config->q);
        }

    } else if (queue_length(q) == 1 && time_until_preempt <= 0 && current_burst_type == CPU_BURST) {
        // 3a) see if any time slices expire (RR mode only) but no other processes are on queue (just print a message)
        current_burst[0] = time_left_in_burst;
        proc->current_burst_start = t;
        printf("time %dms: Time slice expired; no preemption because ready queue is empty [Q <empty>]\n", t);
    } else if (queue_length(q) > 1 && time_until_preempt <= 0 && current_burst_type == CPU_BURST) {
        // 3b) see if any time slices expire (RR mode only); move them to the end of the queue
        current_burst[0] = time_left_in_burst;
        proc->cpu_burst_preempted = true;
        proc->current_burst_type = CX_OFF;
        proc->current_burst_start = t;
        printf("time %dms: Time slice expired; process %c preempted with %dms to go ", t, proc->id, time_left_in_burst);
        print_queue_items(config->q);

    }
}

// 4) see if any processes finish their IO (requeue them)
void RR_handle_io_done(settings *config, process *proc, int t) {
    if (proc->state != BLOCKED || proc->current_burst_type != IO_BURST) return;

    int time_in_burst = t - proc->current_burst_start; 
    int current_burst_length = proc->burst[proc->counter_cpu_burst][1];
    int burst_time_remaining = current_burst_length - time_in_burst;
    if (burst_time_remaining > 0) return;

    // If we get here, the process was blocked due to being in an IO_BURST, which is now done
    proc->state = READY;
    proc->burst[proc->counter_cpu_burst][1] = 0;
    ++proc->counter_cpu_burst;
    proc->current_burst_type = NO_BURST;
    proc->current_burst_start = -1;
    queue_push(config->q, proc, config->rr_queue_push_end);
    printf("time %dms: Process %c completed I/O; added to ready queue ", t, proc->id);
    print_queue_items(config->q);
}

void RR_handle_arrival(settings *config, process *proc, int t) {
    if (proc->state != NOT_YET_ARRIVED) return;
    if (proc->t_arrive > t) return;

    // If we get here, process should arrive now
    proc->state = READY;
    proc->counter_cpu_burst = 0;
    proc->current_burst_type = NO_BURST;
    proc->current_burst_start = t;
    queue_push(config->q, proc, config->rr_queue_push_end);
    printf("time %dms: Process %c arrived; added to ready queue ", t, proc->id);
    print_queue_items(config->q);
}

void RR_clock_tick(settings *config, int t) {
    process** procs = config->procs;
    int num_procs = config->num_procs;
    queue *q = config->q;

    process *current_proc = queue_peek(q);
    // 1a) see if the running process finishes its cx_off switch (either take them off the queue or requeue)
    // 1b) see if the process using the cpu finishes its cx_on switch (mark that the cpu burst has begun)
    // 2) see if the processes using the cpu finishes its burst (either take them off the queue or requeue)
    // 3) see if any time slices expire (RR mode only); move them to the end of the queue
    if (current_proc != NULL) RR_transition_running_process(config, current_proc, t);

    // 4) see if any processes finish their IO (requeue them)
    for (int i = 0; i < num_procs; ++i) RR_handle_io_done(config, procs[i], t);

    // 5) see if any processes arrive (add them to queue after breaking ties)
    for (int i = 0; i < num_procs; ++i) RR_handle_arrival(config, procs[i], t);

    // if the current proc was popped off queue because it CX_OFFed, start the context switch on the next one
    current_proc = queue_peek(q);
    if (current_proc != NULL && current_proc->state == READY) {
        current_proc = queue_peek(q);
        current_proc->state = RUNNING;
        current_proc->current_burst_type = CX_ON;
        current_proc->current_burst_start = t;
        //printf("* time %dms: Process %c found not running at top of queue. Beginning context switch on (will take %d ms). ", t, current_proc->id, config->t_cx / 2);
        //print_queue_items(config->q);
    }

    // 6) see if there are no more processes (exit simulation)
    // No code needed. Checked in for loop
}

void RR(settings *config) {
    //print_config(config);
    for (int i = 0; i < config->num_procs; ++i) print_process_summary(config->procs[i]);
    const char *sim_name = (config->time_slice == INT_MAX) ? "FCFS" : "RR"; 
    printf("time 0ms: Simulator started for %s ", sim_name);
    print_queue_items(config->q);

    int t = 0;
    for (t = 0; num_remaining_procs(config) > 0 || queue_length(config->q) > 0; ++t) {
        RR_clock_tick(config, t);
    }
    printf("time %dms: Simulator ended for %s ", t - 1, sim_name);
    print_queue_items(config->q);

}
