#pragma once

#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h> 

#include "settings.h"
#include "process.h"
#include "queue.h"
void FCFS(settings *config);
void RR(settings *config);
void RR_transition_running_process(settings *config, process *proc, int t);
void RR_handle_io_done(settings *config, process *proc, int t);
void RR_handle_arrival(settings *config, process *proc, int t);
