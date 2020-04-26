#pragma once

#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h> 

#include "settings.h"
#include "process.h"
#include "queue.h"

#ifndef DISPLAY_MAX_T
#define DISPLAY_MAX_T INT_MAX
#endif

void print_event(settings *config, int t, const char *fmt, ...);
void FCFS(settings *config);
void RR(settings *config);
void RR_transition_running_process(settings *config, process *proc, int t);
void RR_handle_io_done(settings *config, process *proc, int t);
void RR_handle_arrival(settings *config, process *proc, int t);
