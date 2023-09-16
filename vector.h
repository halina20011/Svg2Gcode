// Copyright (C) 2023  halina20011
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef VECTOR
#define VECTOR

#include <stdlib.h>

struct Vector{
    double *data;
    size_t size;
    size_t maxSize;
};

#define DEFAULT_MAX_SIZE 20

void newVector(struct Vector **vector){
    *vector = malloc(sizeof(struct Vector));
    (*vector)->maxSize = DEFAULT_MAX_SIZE;
    (*vector)->size = 0;
    (*vector)->data = malloc(sizeof(double) * DEFAULT_MAX_SIZE);
}

void vectorPush(struct Vector *vector, double val){
    if(vector->maxSize <= vector->size){
        vector->maxSize *= 2;
        vector->data = realloc(vector->data, sizeof(double) * vector->maxSize);
    }

    vector->data[vector->size++] = val;
}

bool vectorPop(struct Vector *vector, double *val){
    if(vector->size <= 0){
        return false;
    }
    
    *val = vector->data[vector->size--];
    return true;
}

double *vectorGetDataAndFree(struct Vector *vector, size_t *size){
    double *data = vector->data;
    *size = vector->size / 2;
    free(vector);
    return data;
}

void freeVector(struct Vector *vector){
    free(vector->data);
    free(vector);
}

#endif
