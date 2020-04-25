#include <stdio.h>
#include <string.h>
#include "queue.h"
queue *queue_open() {
    queue *q = malloc(sizeof(queue));
    memset(q, 0, sizeof(queue));
    q->capacity = 26;
    q->items = calloc(sizeof(process*), q->capacity);
    return q;
}

void queue_close(queue *q) {
    free(q->items);
    free(q);
}

void queue_push(queue *q, process *proc, bool to_end) {
    if (to_end) {
        if (q->items[q->start_idx] == NULL) {
            // If queue is empty, place item at start and return
            q->items[q->start_idx] = proc;
            return;
        }

        int empty_idx = (q->start_idx + 1) % q->capacity;
        for (; empty_idx != q->start_idx; empty_idx = (empty_idx + 1) % q->capacity) {
            if (q->items[empty_idx] == NULL) break;
        }
        if (empty_idx == q->start_idx) {
            fprintf(stderr, "Queue got filled somehow.\n");
            exit(EXIT_FAILURE);
        }
        q->items[empty_idx] = proc;
    } else {
        q->start_idx = (q->start_idx + q->capacity - 1) % q->capacity;
        if (q->items[q->start_idx] != NULL) {
            fprintf(stderr, "Queue got filled somehow.\n");
            exit(EXIT_FAILURE);
        }
        q->items[q->start_idx] = proc;
    }
}

process *queue_pop(queue *q) {
    process *result = q->items[q->start_idx];
    q->items[q->start_idx] = NULL;
    q->start_idx = (q->start_idx + 1) % q->capacity;
    return result;
}

process *queue_peek(queue *q) {
    return q->items[q->start_idx];
}

int queue_length(const queue *q) {
    if (q->items[q->start_idx] == NULL) return 0;
    int i = q->start_idx;
    int count = 0;
    do  {
        ++count;
        i = (i + 1) % q->capacity;
    } while (i != q->start_idx && q->items[i] != NULL);
    return count;
}

void print_queue(const queue *q) {
    printf("Processes in Queue (%d, capacity: %d): ", queue_length(q), q->capacity);
    int i = q->start_idx;
    while (q->items[i] != NULL) {
        printf("%c ", q->items[i]->id);
        i = (i + 1) % q->capacity;
    }
    printf("\n");
}
