#include "vec.h"

void f64x3_neg(f64x3_t* a) {
    a->x = -a->x;
    a->y = -a->y;
    a->z = -a->z;
}

void f64x3_add(f64x3_t* a, f64x3_t* b) {
    a->x += b->x;
    a->y += b->y;
    a->z += b->z;
}

void f64x3_sub(f64x3_t* a, f64x3_t* b) {
    a->x -= b->x;
    a->y -= b->y;
    a->z -= b->z;
}

void f64x3_mul(f64x3_t* a, double b) {
    a->x *= b;
    a->y *= b;
    a->z *= b;
}

double f64x3_length(f64x3_t* a) {
    return sqrt(a->x * a->x + a->y * a->y + a->z * a->z);
}

void f64x3_normalize(f64x3_t* a) {
    double rcp_length = 1.0 / sqrt(a->x * a->x + a->y * a->y + a->z * a->z);
    a->x *= rcp_length;
    a->y *= rcp_length;
    a->z *= rcp_length;
}

double f64x3_dot(f64x3_t* a, f64x3_t* b) {
    return a->x * b->x + a->y * b->y + a->z * b->z;
}
