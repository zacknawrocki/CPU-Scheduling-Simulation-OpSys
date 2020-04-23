#include <stdio.h>
#include <string.h>
#include "queue.h"
queue *queue_open() {
    queue *q = malloc(sizeof(queue));
    memset(q, 0, sizeof(queue));
    q->capacity = 26;
    q->items = calloc(sizeof(process*), q->capacity);
}

void queue_close(queue *q) {
    free(q->items);
    free(q);
}

void queue_push(queue *q, process *proc, bool to_end) {
    if (to_end) {
        q->end_idx = (q->end_idx + q->capacity + 1) % q->capacity;
        q->items[q->end_idx] = proc;
    } else {
        q->start_idx = (q->start_idx + q->capacity - 1) % q->capacity;
        q->items[q->start_idx] = proc;
    }
}

process *queue_pop(queue *q) {
    if (q->start_idx == q->end_idx) return NULL;
    process *result = q->items[q->start_idx];
    q->items[q->start_idx] = NULL;
    q->start_idx = (q->start_idx + q->capacity + 1) % q->capacity;
    return result;
}

process *queue_peek(queue *q) {
    if (q->start_idx == q->end_idx) return NULL;
    return q->items[q->start_idx];
}

int queue_length(const queue *q) {
    if (q->end_idx > q->start_idx) {
        return q->end_idx - q->start_idx;
    } else {
        return q->end_idx + (q->capacity - q->start_idx);
    }
}

void print_queue(const queue *q) {
    printf("Processes in Queue (%d): ", queue_length(q));
    bool first = true;
    int i = q->start_idx;
    while (!first && i != q->start_idx) {
        if (!first) ++i;
        printf("%c ", q->items[i]->id);
    }
    printf("\n");
}
