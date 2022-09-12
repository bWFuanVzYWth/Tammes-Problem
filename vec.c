#include "vec.h"

void neg(vec3* a) {
    a->x = -a->x;
    a->y = -a->y;
    a->z = -a->z;
}

void add(vec3* a, vec3* b) {
    a->x += b->x;
    a->y += b->y;
    a->z += b->z;
}

void sub(vec3* a, vec3* b) {
    a->x -= b->x;
    a->y -= b->y;
    a->z -= b->z;
}

void scale(vec3* a, double b) {
    a->x *= b;
    a->y *= b;
    a->z *= b;
}

double length(vec3* a) {
    return sqrt(a->x * a->x + a->y * a->y + a->z * a->z);
}

void normalize(vec3* a) {
    double rcp_length = 1.0 / sqrt(a->x * a->x + a->y * a->y + a->z * a->z);
    a->x *= rcp_length;
    a->y *= rcp_length;
    a->z *= rcp_length;
}

double dot(vec3* a, vec3* b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}
