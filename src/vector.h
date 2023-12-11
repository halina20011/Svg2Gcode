#ifndef VECTOR
#define VECTOR

#include <stdlib.h>
#include <stdbool.h>

#define DEFAULT_MAX_SIZE 20

struct Vector{
    double *data;
    size_t size;
    size_t maxSize;
};

void newVector(struct Vector **vector);
void vectorPush(struct Vector *vector, double val);

bool vectorPop(struct Vector *vector, double *val);
double *vectorGetDataAndFree(struct Vector *vector, size_t *size);
void freeVector(struct Vector *vector);

#endif
