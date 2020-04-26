#pragma once
#include <stdbool.h>
#include "process.h"
#include "queue.h"

typedef struct result {
    double avg_burst_time;
    int total_turnaround_time;
    double avg_turnaround_time;
    int num_cxs;
    int num_preemptions;
} result;

typedef struct settings {
    int num_procs;
    process **procs;
    int t_cx;
    queue *q;
    double alpha;
    int time_slice;
    bool rr_queue_push_end;
    result results;
} settings;


settings *copy_config(const settings *config_tmpl);
void print_config(const settings *config);
void free_config(settings *config);
int num_remaining_procs(const settings *config);
void compute_initial_results(settings *config);
const result *compute_results(settings *config);
