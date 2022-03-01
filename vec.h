#pragma once

#include <math.h>

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

static inline void neg(vec3* a) {
    a->x = -a->x;
    a->y = -a->y;
    a->z = -a->z;
}

static inline void add(vec3* a, vec3* b) {
    a->x += b->x;
    a->y += b->y;
    a->z += b->z;
}

static inline void sub(vec3* a, vec3* b) {
    a->x -= b->x;
    a->y -= b->y;
    a->z -= b->z;
}

static inline void scale(vec3* a, double b) {
    a->x *= b;
    a->y *= b;
    a->z *= b;
}

static inline void normalize(vec3* a) {
    double rcp_length = 1.0 / sqrt(a->x * a->x + a->y * a->y + a->z * a->z);
    a->x *= rcp_length;
    a->y *= rcp_length;
    a->z *= rcp_length;
}

static inline double dot(vec3* a, vec3* b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}
