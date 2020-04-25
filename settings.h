#pragma once
#include <stdbool.h>
#include "process.h"
#include "queue.h"

typedef struct settings {
    int num_procs;
    process **procs;
    int t_cx;
    queue *q;
    //int preemptions;
    double alpha;
    int time_slice;
    bool rr_queue_push_end;
} settings;


settings *copy_config(const settings *config_tmpl);
void print_config(const settings *config);
void free_config(settings *config);
int num_remaining_procs(const settings *config);
