#include <math.h>
#include <stdint.h>
#include <time.h>

#include <stdio.h>

#include "random_point.h"

// PCG伪随机数生成器，因为<stdlib.h>的rand()不支持同时使用多个种子

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

uint32_t pcg32_random_r(pcg32_random_t* rng) {
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// -1.0 <= random_1 < 1.0
double random_1(pcg32_random_t* rng) {
    return pcg32_random_r(rng) * (1.0 / 0x80000000) - 1.0;
}

// 在单位圆内随机生成一个点
void circle_point_picking(double* x, double* y, double* r2, pcg32_random_t* rng) {
    do {
        *x = random_1(rng);
        *y = random_1(rng);
        *r2 = (*x) * (*x) + (*y) * (*y);
    } while (*r2 >= 1.0);
}

// 在单位球面上随机生成一个点
void sphere_point_picking(f64x3_t* v, pcg32_random_t* rng) {
    double x, y, r2, tmp;
    circle_point_picking(&x, &y, &r2, rng);
    tmp = 2.0 * sqrt(1.0 - r2);
    v->x = x * tmp;
    v->y = y * tmp;
    v->z = -2.0 * r2 + 1.0;
}
