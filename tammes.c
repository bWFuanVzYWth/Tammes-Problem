#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mt19937ar.h"

typedef struct {
    double x;
    double y;
    double z;
} vec3;

vec3 add(vec3 a, vec3 b) {
    vec3 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}

vec3 mul(double a, vec3 b) {
    vec3 c;
    c.x = a * b.x;
    c.y = a * b.y;
    c.z = a * b.z;
    return c;
}

double dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3 normalize(vec3 a) {
    double rcp_length = 1.0 / sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    vec3 b;
    b.x = a.x * rcp_length;
    b.y = a.y * rcp_length;
    b.z = a.z * rcp_length;
    return b;
}

#define POINT_LIST_LENGTH 10000  // Enough points. Fejes Tóth, 1943

vec3 sphere_point_picking(void) {
    // https://mathworld.wolfram.com/SpherePointPicking.html
    double x1, x2, tmp_1, tmp_2;
    while (1) {
        x1 = genrand_res53() * 2.0 - 1.0;  // (-1,1)
        x2 = genrand_res53() * 2.0 - 1.0;  // (-1,1)
        tmp_1 = x1 * x1 + x2 * x2;
        if (tmp_1 < 1.0)
            break;
    }

    vec3 a;
    // a.x = 2.0 * x1 * sqrt(1.0 - x1 * x1 - x2 * x2);
    // a.y = 2.0 * x2 * sqrt(1.0 - x1 * x1 - x2 * x2);
    // a.z = 1.0 - 2.0 * (x1 * x1 + x2 * x2);
    tmp_2 = 2.0 * sqrt(1.0 - tmp_1);
    a.x = x1 * tmp_2;
    a.y = x2 * tmp_2;
    a.z = 1.0 - 2.0 * tmp_1;

    return a;
}

void random_point(vec3 point[POINT_LIST_LENGTH], size_t point_num) {
    for (size_t i = 0; i < point_num; i++) {
        point[i] = sphere_point_picking();
    }
}

void output(vec3 point[POINT_LIST_LENGTH], size_t point_num) {
    putchar('{');
    for (size_t i = 0; i < point_num; i++) {
        printf("{%1.6lf,%1.6lf,%1.6lf}", point[i].x, point[i].y, point[i].z);
        if (i < point_num - 1)
            putchar(',');
    }
    putchar('}');
}

void copy(vec3 best_point[POINT_LIST_LENGTH], vec3 point[POINT_LIST_LENGTH], size_t point_num) {
    for (size_t i = 0; i < point_num; i++) {
        best_point[i] = point[i];
    }
}

double find_max_dot(vec3 point[POINT_LIST_LENGTH], size_t point_num) {
    double max_dot = -2.0;
    for (size_t i = 0; i < point_num - 1; i++) {
        for (size_t j = i + 1; j < point_num; j++) {
            double tmp_dot = dot(point[i], point[j]);
            if (tmp_dot > max_dot) {
                max_dot = tmp_dot;
            }
        }
    }
    return max_dot;
}

void find_nearest_point(vec3 point[POINT_LIST_LENGTH], size_t point_num, size_t* a, size_t* b) {
    double max_dot = -2.0;
    for (size_t i = 0; i < point_num - 1; i++) {
        for (size_t j = i + 1; j < point_num; j++) {
            double tmp_dot = dot(point[i], point[j]);
            if (tmp_dot > max_dot) {
                max_dot = tmp_dot;
                *a = i;
                *b = j;
            }
        }
    }
}

vec3 get_move_vec(double temperature) {
    vec3 move_dir = sphere_point_picking();
    double move_len = temperature * sqrt(log(1.0 / (1.0 - genrand_res53())));  // TODO T
    return mul(move_len, move_dir);
}

void annealing(vec3 point[POINT_LIST_LENGTH], size_t point_num) {
    const size_t step = 100000000;
    const double speed = 10.0;

    double best_max_dot = 2.0;

    double temperature = 1.0;

    for (size_t i = 0; i < step; i++) {
        temperature = exp(i * (-speed / step));

        size_t a, b;
        find_nearest_point(point, point_num, &a, &b);
        best_max_dot = dot(point[a], point[b]);
        vec3 new_point_a = normalize(add(point[a], get_move_vec(temperature)));
        vec3 new_point_b = normalize(add(point[b], get_move_vec(temperature)));

        if (best_max_dot > dot(new_point_a, new_point_b) || genrand_res53() < temperature) {
            point[a] = new_point_a;
            point[b] = new_point_b;
        }
    }
    printf("min_angle = %lf\n", 180.0 * acos(best_max_dot) / M_PI);
}

int main(void) {
    vec3 point[POINT_LIST_LENGTH];
    size_t point_num = 60;

    init_genrand((unsigned int)time(NULL));

    random_point(point, point_num);
    annealing(point, point_num);
    output(point, point_num);

    return 0;
}