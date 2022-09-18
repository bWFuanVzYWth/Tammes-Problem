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

int check_input(uint32_t* point_num, uint32_t* mode) {
    if (*point_num < 4) {
        fprintf(stderr, "节点数不能小于4\n");
        return -1;
    }
    switch (*mode) {
        case 1:
            if (*point_num % 1 == 0)
                break;
        case 2:
            if (*point_num % 2 == 0)
                break;
        default:
            *mode = 1;
            fprintf(stderr, "参数错误，已忽略假设，并将此参数设置为默认值1\n假设对称性暂时只支持：1不对称 2中心对称\n其中中心对称要求节点数为2的倍数才能启用\n");
            break;
    }
    return 0;
}

void output_mma(object_t* object) {
    char filename[127] = {0};
    sprintf(filename, "%u_%lf_mma.txt", object->point_num, object->angle);
    FILE* fp = fopen(filename, "wb");
    fprintf(fp, "point={");
    for (int i = 0; i < object->point_num; i++) {
        fprintf(fp, "{%1.16lf,%1.16lf,%1.16lf}%c", object->pos[i].x, object->pos[i].y, object->pos[i].z, (i == object->point_num - 1) ? '}' : ',');
    }
    fclose(fp);
    free(object->pos);
}

void output_Aqi(object_t* object) {
    char filename[127] = {0};
    sprintf(filename, "%u_%lf_Aqi.txt", object->point_num, object->angle);
    FILE* fp = fopen(filename, "wb");
    fprintf(fp, "AqiDYBP=");
    for (int i = 0; i < object->point_num; i++) {
        fprintf(fp, "%1.16lf,%1.16lf,%1.16lf%c", object->pos[i].x, object->pos[i].y, object->pos[i].z, (i == object->point_num - 1) ? '|' : '\0');
    }
    fclose(fp);
    free(object->pos);
}

void output_csv(object_t* object) {
    char filename[127] = {0};
    sprintf(filename, "%u_%lf.csv", object->point_num, object->angle);
    FILE* fp = fopen(filename, "wb");
    for (int i = 0; i < object->point_num; i++) {
        fprintf(fp, "%1.16lf,%1.16lf,%1.16lf\n", object->pos[i].x, object->pos[i].y, object->pos[i].z);
    }
    fclose(fp);
    free(object->pos);
}

int main(void) {
    // const i32x3_t version = {0, 2, 5};
    time_t seed_0 = time(NULL);

    uint32_t point_num;
    uint64_t iteration;
    uint32_t repeat;
    uint32_t mode;

    fprintf(stderr, "请依次输入节点数，迭代次数，重试次数，假设对称性，用空格分隔，然后按回车（示例：130 1000000 128 1）：\n");
    scanf("%" PRIu32 " %" PRIu64 " %" PRIu32 " %" PRIu32, &point_num, &iteration, &repeat, &mode);
    if (check_input(&point_num, &mode))
        return -1;

    double time_0 = omp_get_wtime();

    object_t* object = calloc(repeat, sizeof(object_t));
    for (int i = 0; i < repeat; i++) {
        object_init(&object[i], iteration, point_num);
    }

    for (int i = 0; i < repeat; i++) {
        pcg32_random_t rng = {0, seed_0 + i};
        switch (mode) {
            case 2:
                for (int j = 0; j < object[i].point_num / 2; j++) {
                    sphere_point_picking(&(object[i].pos[2 * j]), &rng);
                    f64x3_mov(&(object[i].pos[2 * j + 1]), &(object[i].pos[2 * j]));
                    f64x3_neg(&(object[i].pos[2 * j + 1]));
                }
                break;
            default:
                for (int j = 0; j < object[i].point_num; j++) {
                    sphere_point_picking(&(object[i].pos[j]), &rng);
                }
                break;
        }
    }

#pragma omp parallel for
    for (int i = 0; i < repeat; i++) {
        object[i].angle = tammes(object[i].pos, object[i].point_num, object[i].iteration, mode);
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
    output_csv(&object[best_index]);

    for (int i = 0; i < repeat; i++) {
        object_free(&object[i]);
    }

    return 0;
}