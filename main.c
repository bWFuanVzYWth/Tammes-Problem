#include <inttypes.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "random_point.h"
#include "tammes.h"
#include "vec.h"

#define OUTPUT_PRECISION "1.16"

typedef struct {
    uint64_t iteration;
    uint32_t point_num;
    double angle;
    f64x3_t* pos;
} object_t;

void object_init(object_t* object, uint64_t iteration, uint32_t point_num) {
    object->iteration = iteration;
    object->point_num = point_num;
    object->pos = calloc(point_num, sizeof(f64x3_t));
}

void object_free(object_t* object) {
    free(object->pos);
}

int main(void) {
    // const i32x3_t version = {0, 2, 5};
    time_t seed_0 = time(NULL);

    uint32_t point_num;
    uint64_t iteration;
    uint32_t repeat;

    fprintf(stderr, "请依次输入节点数，迭代次数，重试次数，用空格分隔，然后按回车（示例：130 1000000 128）：\n");
    scanf("%" PRIu32 " %" PRIu64 " %" PRIu32, &point_num, &iteration, &repeat);

    double time_0 = omp_get_wtime();

    object_t* object = calloc(repeat, sizeof(object_t));
    for (int i = 0; i < repeat; i++) {
        object_init(&object[i], iteration, point_num);
        printf("id=%d, iteration=%" PRIu64 ", point_num=%" PRIu32 "\n", i, object[i].iteration, object[i].point_num);
    }

    for (int i = 0; i < repeat; i++) {
        pcg32_random_t rng = {0, seed_0 + i};
        for (int j = 0; j < object[i].point_num; j++) {
            sphere_point_picking(&(object[i].pos[j]), &rng);
        }
    }

#pragma omp parallel for
    for (int i = 0; i < repeat; i++) {
        object[i].angle = tammes(object[i].pos, object[i].point_num, object[i].iteration);
        fprintf(stderr, "线程%4d优化完毕，最小夹角=%1.6lf\n", i, object[i].angle);
    }

    double adv = 0.0;
    uint32_t best_index = 0;
    for (int i = 0; i < repeat; i++) {
        adv += object[i].angle;
        if (object[i].angle > object[best_index].angle)
            best_index = i;
    }
    adv /= repeat;

    double time_1 = omp_get_wtime();
    double time = time_1 - time_0;
    fprintf(stderr, "所有%" PRIu32 "个线程优化完毕，平均%lf，用时%lfs，被最大化的最小夹角=%1.16lf\n", repeat, adv, time, object[best_index].angle);

    // 将坐标输出为蓝图

    for (int i = 0; i < repeat; i++) {
        object_free(&object[i]);
    }

    return 0;
}