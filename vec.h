#pragma once

#include <math.h>
#include <stdint.h>

typedef struct {
    int x;
    int y;
    int z;
} int3;

typedef struct {
    double x;
    double y;
    double z;
} vec3;

void neg(vec3* a);

void add(vec3* a, vec3* b);
void sub(vec3* a, vec3* b);
void scale(vec3* a, double b);

double length(vec3* a);
double dot(vec3* a, vec3* b);

void normalize(vec3* a);
