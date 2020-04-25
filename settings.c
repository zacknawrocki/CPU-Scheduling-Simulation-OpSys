#include <string.h>
#include <stdio.h>
#include <limits.h>
#include "settings.h"
#include "process.h"
settings *copy_config(const settings *config_tmpl) {
    settings *config = malloc(sizeof(settings));
    memcpy(config, config_tmpl, sizeof(settings));

    config->procs = calloc(config->num_procs, sizeof(process));
    for (int i = 0; i < config->num_procs; ++i) config->procs[i] = copy_process(config_tmpl->procs[i]);
    config->q = queue_open(); 
    return config;
}

void free_config(settings *config) {
    for (int i = 0; i < config->num_procs; ++i) free_process(config->procs[i]);
    free(config->procs);
    queue_close(config->q);
    free(config);
}


void print_config(const settings *config) {
    printf("Number of Processes: %d\n", config->num_procs);
    printf("Processes: \n");
    for (int i = 0; i < config->num_procs; ++i) print_process(config->procs[i]);
    printf("Context Switch Time: %d (Half: %d)\n", config->t_cx, config->t_cx/2);
    print_queue(config->q);

    printf("Alpha: %f\n", config->alpha);
    if (config->time_slice == INT_MAX) {
        printf("Time Slice: infinite\n");;
    } else {
        printf("Time Slice: %d\n", config->time_slice);
    }
    printf("Round Robbin Queue Pushes to End: %s\n", config->rr_queue_push_end ? "yes" : "no");
}
