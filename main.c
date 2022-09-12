#define _USE_MATH_DEFINES
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "vec.h"
#include "dSFMT.h"

#include "tammes.h"

// 在单位圆内随机生成一个点
void circle_point_picking(double* x1, double* x2, double* r2, dsfmt_t* dsfmt) {
    do {
        *x1 = dsfmt_genrand_open_open(dsfmt) * 2.0 - 1.0;
        *x2 = dsfmt_genrand_open_open(dsfmt) * 2.0 - 1.0;
        *r2 = (*x1) * (*x1) + (*x2) * (*x2);
    } while (*r2 >= 1.0);
}

// 在单位球面上随机生成一个点
void sphere_point_picking(vec3* v, dsfmt_t* dsfmt) {
    double x1, x2, r2, tmp;
    circle_point_picking(&x1, &x2, &r2, dsfmt);
    tmp = 2.0 * sqrt(1.0 - r2);
    v->x = x1 * tmp;
    v->y = x2 * tmp;
    v->z = -2.0 * r2 + 1.0;
}

void creat_random_object(object* object, uint32_t seed, uint32_t point_num) {
    dsfmt_t dsfmt;
    object->point_num = point_num;
    object->point = calloc(point_num, sizeof(list));
    object->best_point = calloc(point_num, sizeof(list));
    dsfmt_init_gen_rand(&dsfmt, seed);
    for (uint32_t i = 0; i < point_num; i++) {
        sphere_point_picking(&object->point[i].pos, &dsfmt);
        object->point[i].id = i;
        object->point[i].next = NULL;
    }
}

#define OUTPUT_PRECISION "1.16"

// debug，其中的坐标是mathematica中的列表格式
void output_to_debug(int3 version, time_t seed_0, object* object_list, uint32_t best_index, uint32_t point_num, uint64_t iteration, uint32_t repeat, double time, double adv) {
    char filename[256] = {0};
    sprintf(filename, "%u_%1.6lf_debug.txt", point_num, object_list[best_index].angle);
    FILE* fp = fopen(filename, "w");

    // 出错时重定向输出到屏幕
    if (fp == NULL)
        fp = stderr;

    fprintf(fp, "{");
    for (uint32_t i = 0; i < object_list[best_index].point_num; i++) {
        if (i < object_list[best_index].point_num - 1)
            fprintf(fp, "{%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf},", object_list[best_index].best_point[i].pos.x, object_list[best_index].best_point[i].pos.y, object_list[best_index].best_point[i].pos.z);
        else
            fprintf(fp, "{%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf}}\n", object_list[best_index].best_point[i].pos.x, object_list[best_index].best_point[i].pos.y, object_list[best_index].best_point[i].pos.z);
    }

    fprintf(fp, "version = %d.%d.%d\n", version.x, version.y, version.z);
    fprintf(fp, "seed_0 = %llu, point_num = %u, iteration = %llu, repeat = %u\n", seed_0, point_num, iteration, repeat);
    fprintf(fp, "best_index = %u, best_angle = %" OUTPUT_PRECISION "lf, time = %lf, adv = %lf\n", best_index, object_list[best_index].angle, time, adv);
    for (int i = 0; i < repeat; i++)
        fprintf(fp, "id = %4d , angle = %" OUTPUT_PRECISION "lf\n", i, object_list[i].angle);

    if (fp == stderr) {
        for (int i = 0; i < 3; i++)
            fprintf(stderr, "!!!写入文件时发生错误，尝试将结果输出至屏幕，请手动保存后再关闭程序!!!\n");
        while (1)
            getchar();
    } else {
        fclose(fp);
    }
}

// 输出为可使用mod直接导入的蓝图格式
void output_to_AqiDYBP(object* object_list, uint32_t best_index, uint32_t point_num) {
    char filename[256] = {0};
    sprintf(filename, "%u_%1.6lf.txt", point_num, object_list[best_index].angle);
    FILE* fp = fopen(filename, "w");

    // 出错时直接返回
    if (fp == NULL)
        return;

    fprintf(fp, "AqiDYBP=");
    for (uint32_t i = 0; i < object_list[best_index].point_num; i++) {
        if (i < object_list[best_index].point_num - 1)
            fprintf(fp, "%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf|", object_list[best_index].best_point[i].pos.x, object_list[best_index].best_point[i].pos.y, object_list[best_index].best_point[i].pos.z);
        else
            fprintf(fp, "%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf", object_list[best_index].best_point[i].pos.x, object_list[best_index].best_point[i].pos.y, object_list[best_index].best_point[i].pos.z);
    }

    fclose(fp);
}

int main(void) {
    const int3 version = {0, 2, 5};
    time_t seed_0 = time(NULL);

    uint32_t point_num;
    uint64_t iteration;
    uint32_t repeat;

    fprintf(stderr, "请依次输入节点数，迭代次数，重试次数，用空格分隔，然后按回车（示例：130 1000000 128）：\n");
    scanf("%" PRIu32 " %" PRIu64 " %" PRIu32, &point_num, &iteration, &repeat);

    double time_0 = omp_get_wtime();

    object* object_list = calloc(repeat, sizeof(object));

    // 生产随机初始状态
    for (int i = 0; i < repeat; i++) {
        creat_random_object(object_list + i, (uint32_t)seed_0 + i, point_num);
    }

    // 对每个初始状态分别开一个线程进行优化
#pragma omp parallel for
    for (int i = 0; i < repeat; i++) {
        tammes(object_list + i, iteration);
        fprintf(stderr, "线程%4d优化完毕，最小夹角=%1.6lf\n", i, object_list[i].angle);
    }

    // 找出最好的那一个，并评估收敛水平
    double adv = 0.0;
    uint32_t best_index = 0;
    for (int i = 0; i < repeat; i++) {
        adv += object_list[i].angle;
        if (object_list[i].angle > object_list[best_index].angle)
            best_index = i;
    }
    adv /= repeat;

    double time_1 = omp_get_wtime();
    double time = time_1 - time_0;
    fprintf(stderr, "所有%" PRIu32 "个线程优化完毕，平均%lf，用时%lfs，被最大化的最小夹角=%1.16lf\n", repeat, adv, time, object_list[best_index].angle);

    // 将坐标输出为蓝图
    output_to_AqiDYBP(object_list, best_index, point_num);
    output_to_debug(version, seed_0, object_list, best_index, point_num, iteration, repeat, time, adv);

    for (int i = 0; i < repeat; i++) {
        free(object_list[i].point);
        free(object_list[i].best_point);
    }
    free(object_list);

    return 0;
}