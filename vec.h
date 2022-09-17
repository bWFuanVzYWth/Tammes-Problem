#ifndef _VEC_HEAD_
#define _VEC_HEAD_

#include <math.h>
#include <stdint.h>

typedef struct f64x3_t f64x3_t;
typedef struct i32x3_t i32x3_t;

struct i32x3_t {
    int32_t x;
    int32_t y;
    int32_t z;
};

struct f64x3_t {
    double x;
    double y;
    double z;
};

void f64x3_neg(f64x3_t* a);

void f64x3_add(f64x3_t* a, f64x3_t* b);
void f64x3_sub(f64x3_t* a, f64x3_t* b);
void f64x3_mul(f64x3_t* a, double b);

double f64x3_length(f64x3_t* a);
double f64x3_dot(f64x3_t* a, f64x3_t* b);

void f64x3_normalize(f64x3_t* a);

#endif
