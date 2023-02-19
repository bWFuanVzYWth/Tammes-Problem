#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "random_point.h"

typedef enum error_type {
    no_error = 0,
    no_argument,
    unknow_argument,
    wrong_formart,
    illegal_point_num,
    illegal_symmetry,
    file_no_found,
    out_of_memory,
} error_t;

#ifdef ENGLISH_DOC

#define HELP_DOC \
    "https://github.com/bWFuanVzYWth/Tammes-Problem v3.0 MIT license\n\
Usage: tammes.exe [arg1] [arg2] ...\n\
Arguments:\n\
-N=<num>    Number of points, n>=2. Required if no -F.\n\
-F=<path>   Initialization points by reading from file instead of random.\n\
-M=<1,2>    Assume that the points has m-fold symmetry. default: 1.\n\
-S=<num>    Specify random seed. default: current timestamp, in ns.\n\
Examples:\n\
tammes.exe -N=130\n\
tammes.exe -N=130 -M=2 -S=114514\n\
tammes.exe -F=C:\\tammes\\input.csv"

char* err_tex[] = {
    "no error",           // no_error
    "no argument",        // no_argument
    "unknow argument",    // unknow_argument
    "wrong formart",      // wrong_formart
    "illegal point num",  // illegal_point_num
    "illegal symmetry",   // illegal_symmetry
    "file no found",      // file_no_found
    "out of memory"       // out_of_memory
};

#define ERROR_DOC "!!! Unable to start: %s !!!\n"

#else  // Chinese

#define HELP_DOC \
    "https://github.com/bWFuanVzYWth/Tammes-Problem v3.0 MIT license\n\
使用方式: tammes.exe [参数1] [参数2] ...\n\
参数:\n\
-N=<num>    球面上点的个数, N>=2. 必填，除非启用-F选项\n\
-F=<path>   从文件读取初始坐标，而不是从随机状态开始\n\
-M=<1,2>    假定点的排列具有M重对称性. 只能填1或2, 默认值=1\n\
-S=<num>    指定伪随机种子. 默认值=精确到纳秒的当前时间.\n\
示例:\n\
tammes.exe -N=130\n\
tammes.exe -N=130 -M=2 -S=114514\n\
tammes.exe -F=C:\\tammes\\input.csv"

char* err_tex[] = {
    "没有发现异常",        // no_error
    "缺少启动参数",        // no_argument
    "未知启动参数",        // unknow_argument
    "错误的参数格式",      // wrong_formart
    "不合理的点数",        // illegal_point_num
    "不合理的假定对称性",  // illegal_symmetry
    "找不到文件",          // file_no_found
    "无法分配内存"         // out_of_memory
};

#define ERROR_DOC "!!! 无法启动: %s !!!\n"

#endif

#define panic(err)                                    \
    {                                                 \
        if (err) {                                    \
            fprintf(stderr, ERROR_DOC, err_tex[err]); \
            error_code = err;                         \
            goto need_help;                           \
        }                                             \
    }

#define MAX_PATH 32768

typedef enum bool_type { false = 0, true = 1 } bool;

typedef struct config_type {
    char file_path[MAX_PATH];
    uint64_t seed;
    int64_t point_num;
    int64_t symmetry;
    bool use_seed;
    bool read_from_file;
} config_t;

uint64_t get_timestamp(void) {
    struct timespec t;
    clock_gettime(0, &t);
    return (uint64_t)t.tv_sec * 1000000000 + (uint64_t)t.tv_nsec;
}

error_t init_from_file(FILE* fp, f64x3_t** point, int64_t* point_num) {
    (*point_num) = 0;
    f64x3_t tmp;
    while (fscanf(fp, "%lf,%lf,%lf", &tmp.x, &tmp.y, &tmp.z) == 3) {
        (*point_num)++;
        (*point) = (f64x3_t*)realloc((*point), sizeof(f64x3_t) * (*point_num));
        if ((*point) == NULL)
            return out_of_memory;
        memcpy(&(*point)[(*point_num) - 1], &tmp, sizeof(f64x3_t));
    }
    return no_error;
}

error_t init_from_random(pcg32_random_t* pcg, f64x3_t** point, int64_t point_num) {
    (*point) = (f64x3_t*)realloc((*point), sizeof(f64x3_t) * point_num);
    if ((*point) == NULL)
        return out_of_memory;
    for (int i = 0; i < point_num; i++)
        sphere_point_picking(&(*point)[i], pcg);
    return no_error;
}

void dump_to_csv(FILE* fpw, f64x3_t* point, int64_t point_num) {
    for (int i = 0; i < point_num; i++) {
        fprintf(fpw, "%1.18lf,%1.18lf,%1.18lf\n", point[i].x, point[i].y, point[i].z);
    }
}

error_t check_no_argument(int argc) {
    return argc >= 2 ? no_error : no_argument;
}

error_t check_illegal_point_num(int64_t point_num) {
    return point_num >= 2 ? no_error : illegal_point_num;
}

error_t check_illegal_symmetry(int64_t symmetry) {
    return (symmetry == 1 || symmetry == 2) ? no_error : illegal_symmetry;
}

error_t check_file_no_found(FILE* fp) {
    return fp != NULL ? no_error : file_no_found;
}

error_t check_wrong_formart(int n) {
    return n == 1 ? no_error : wrong_formart;
}

error_t parsing_format(char* argv, config_t* cfg) {
    char argument = 0;
    char string[MAX_PATH] = { 0 };
    sscanf(argv, "-%c=%s", &argument, string);
    switch (argument) {
    case 'N':
        if (cfg->read_from_file == false)
            return check_wrong_formart(sscanf(string, "%" PRId64, &cfg->point_num));
    case 'F':
        cfg->read_from_file = true;
        strcpy(cfg->file_path, string);
        return no_error;
    case 'M':
        return check_wrong_formart(
            sscanf(string, "%" PRId64, &cfg->symmetry));
    case 'S':
        cfg->use_seed = true;
        return check_wrong_formart(sscanf(string, "%" PRIu64, &cfg->seed));
    default:
        return unknow_argument;
    }
}

int main(int argc, char* argv[]) {
    config_t cfg = { {0}, 0, -1, 1, false, false };
    error_t error_code = no_error;
    panic(check_no_argument(argc));
    panic(check_illegal_symmetry(cfg.symmetry));

    for (int i = 1; i < argc; i++)
        panic(parsing_format(argv[i], &cfg));

    f64x3_t* point = NULL;

    if (cfg.read_from_file == true) {
        FILE* fpr = fopen(cfg.file_path, "r");
        panic(check_file_no_found(fpr));
        panic(init_from_file(fpr, &point, &(cfg.point_num)));
        panic(check_illegal_point_num(cfg.point_num));
        fclose(fpr);
    }
    else {
        panic(check_illegal_point_num(cfg.point_num));
        pcg32_random_t pcg = { cfg.use_seed ? cfg.seed : get_timestamp(), 0 };
        panic(init_from_random(&pcg, &point, cfg.point_num));
    }

    // tammes();

    FILE* fpw = fopen("tmp.csv", "wb");
    dump_to_csv(fpw, point, cfg.point_num);
    fclose(fpw);

    return 0;

need_help:
    puts(HELP_DOC);
    return error_code;
}