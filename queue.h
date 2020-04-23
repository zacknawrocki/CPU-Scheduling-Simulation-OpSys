#pragma once
#include <stdlib.h>
#include <stdbool.h>
#include "process.h"

typedef struct queue {
    int capacity;
    int start_idx;
    int end_idx;
    process **items;
} queue;

int queue_length(const queue *q);
process *queue_peek(queue *q);
process *queue_pop(queue *q);
void queue_push(queue *q, process *proc, bool to_end);
void queue_close(queue *q);
queue *queue_open();
void print_queue(const queue *q);
