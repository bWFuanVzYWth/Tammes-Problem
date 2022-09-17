#ifndef _RANDOM_POINT_HEAD_
#define _RANDOM_POINT_HEAD_

#include "vec.h"

typedef struct {
    uint64_t state;
    uint64_t inc;
} pcg32_random_t;

void sphere_point_picking(f64x3_t* v, pcg32_random_t* rng);

#endif
