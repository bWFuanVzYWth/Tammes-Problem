#define _USE_MATH_DEFINES
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "avlmini.h"
#include "dSFMT.h"
#include "vec.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

// struct list* next;     // 8   链表指针
// vec3 pos;              // 8*3 点的坐标
// size_t hash;           // 8   散列表
// int3 hash3_min;        // 4*3 空间划分网格中最小值
// int3 hash3_min;        // 4*3 空间划分网格中最大值
// size_t id;             // 4   点的编号
typedef struct list list;

struct list {
    list* next;
    vec3 pos;
    size_t hash;
    int3 hash3_min;
    int3 hash3_max;
    uint32_t id;
};

// list* point;         // 32*n   点的坐标列表
// double angle;        // 8      最小点对的夹角
// dsfmt_t dsfmt;       // 16*n+4 伪随机发生器状态
// uint32_t point_num;  // 4      点的数量
typedef struct {
    list* point;
    double angle;
    dsfmt_t dsfmt;
    uint32_t point_num;
} object;

// list* point1;          // 8     序号更小的点指针
// list* point2;          // 8     序号更大的点指针
// double cos_angle;      // 8     两点之间的夹角余弦
// struct avl_node node;  // 8*3+4 AVL树的指针
typedef struct {
    list* point1;
    list* point2;
    double cos_angle;
    struct avl_node node;
} tree;

// 对于球面上的n点，总存在两个点，其距离<=D (Fejes Tóth, 1943)
double get_D(uint32_t point_num) {
    double tmp = 1.0 / sin((M_PI * point_num) / (6 * (point_num - 2)));
    return sqrt(4.0 - tmp * tmp);
}

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

size_t to_hash(int x, int y, int z, const uint32_t N) {
    return (x * N + y) * N + z;
}

void refresh_hash(list* point, uint32_t N) {
    int3 hash3;
    // 将空间坐标映射到空间划分网格坐标
    hash3.x = (int)((point->pos.x * 0.5 + 0.5) * N);
    hash3.y = (int)((point->pos.y * 0.5 + 0.5) * N);
    hash3.z = (int)((point->pos.z * 0.5 + 0.5) * N);
    // 计算空间划分网格的邻域范围
    point->hash3_min.x = max(hash3.x - 1, 0);
    point->hash3_min.y = max(hash3.y - 1, 0);
    point->hash3_min.z = max(hash3.z - 1, 0);
    point->hash3_max.x = min(hash3.x + 1, N - 1);
    point->hash3_max.y = min(hash3.y + 1, N - 1);
    point->hash3_max.z = min(hash3.z + 1, N - 1);
    // 将空间划分网格坐标映射到hash
    point->hash = to_hash(hash3.x, hash3.y, hash3.z, N);
}

// AVL树的比较函数，所有节点存放都在一个数组里，地址固定，因此地址的大小作为备选的比较方式
int tree_compare(const void* v_p1, const void* v_p2) {
    tree* p1 = (tree*)v_p1;
    tree* p2 = (tree*)v_p2;

    // 先检查是否为同一节点，避免浮点误差
    if (p1 == p2)
        return 0;
    // 距离不同时的比较，使用math.h中的函数比较浮点数大小，cos越大夹角越小所以是反的
    if (isless(p2->cos_angle, p1->cos_angle))
        return -1;
    if (isless(p1->cos_angle, p2->cos_angle))
        return 1;
    // 相等时的处理，应该很少走到这里
    if (p1 < p2)
        return -1;
    else
        return 1;
}

void point_add(list* point, list** hashmap, const uint32_t N) {
    point->next = hashmap[point->hash];
    hashmap[point->hash] = point;
}

void point_remove(list* point, list** hashmap, const uint32_t N) {
    list* last = (list*)&hashmap[point->hash];
    list* here = hashmap[point->hash];
    while (point != here) {
        last = here;
        here = here->next;
    }
    last->next = here->next;
}

void distance_remove(list* point, list** hashmap, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D) {
    for (size_t hash_x = point->hash3_min.x; hash_x <= point->hash3_max.x; hash_x++) {
        for (size_t hash_y = point->hash3_min.y; hash_y <= point->hash3_max.y; hash_y++) {
            for (size_t hash_z = point->hash3_min.z; hash_z <= point->hash3_max.z; hash_z++) {
                list* p_del = hashmap[to_hash(hash_x, hash_y, hash_z, N)];
                while (p_del != NULL) {
                    double cos_angle = dot(&point->pos, &p_del->pos);
                    if (isgreaterequal(cos_angle, cos_D)) {
                        list* point1 = min(point, p_del);
                        list* point2 = max(point, p_del);
                        tree* p = &distance_list[point1->id * point_num + point2->id];
                        avl_tree_remove(distance, p);
                    }
                    p_del = p_del->next;
                }
            }
        }
    }
}

void distance_add(list* point, list** hashmap, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D) {
    for (size_t hash_x = point->hash3_min.x; hash_x <= point->hash3_max.x; hash_x++) {
        for (size_t hash_y = point->hash3_min.y; hash_y <= point->hash3_max.y; hash_y++) {
            for (size_t hash_z = point->hash3_min.z; hash_z <= point->hash3_max.z; hash_z++) {
                list* p_add = hashmap[to_hash(hash_x, hash_y, hash_z, N)];
                while (p_add != NULL) {
                    double cos_angle = dot(&point->pos, &p_add->pos);
                    if (isgreaterequal(cos_angle, cos_D)) {
                        list* point1 = min(point, p_add);
                        list* point2 = max(point, p_add);
                        tree* p = &distance_list[point1->id * point_num + point2->id];
                        p->point1 = point1;
                        p->point2 = point2;
                        p->cos_angle = cos_angle;
                        avl_tree_add(distance, p);
                    }
                    p_add = p_add->next;
                }
            }
        }
    }
}

void move_point(list* point, vec3* move_vec, list** hashmap, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D) {
    // 从网格中移除这个点，然后从树中移除相关的点对
    point_remove(point, hashmap, N);
    distance_remove(point, hashmap, distance, distance_list, N, point_num, cos_D);
    // 移动点
    add(&point->pos, move_vec);
    normalize(&point->pos);
    // 更新坐标，把相关的点对加入树，然后把这个点加入网格
    refresh_hash(point, N);
    distance_add(point, hashmap, distance, distance_list, N, point_num, cos_D);
    point_add(point, hashmap, N);
}

void tammes(object* object, uint64_t iteration) {
    list* point = object->point;
    uint32_t point_num = object->point_num;

    const double D = get_D(point_num);
    const uint32_t N = (uint32_t)floor(2.0 / D);
    const double cos_D = 1.0 - 0.5 * D * D - 1e-7;

    // 初始化AVL树，然后创建距离矩阵，同时也是AVL树的节点
    struct avl_tree distance;
    avl_tree_init(&distance, &tree_compare, sizeof(tree), AVL_OFFSET(tree, node));
    tree* distance_list = (tree*)calloc(point_num * point_num, sizeof(tree));
    // 创建空间均匀划分网格，可以理解成hashmap
    list** hashmap = (list**)calloc(N * N * N, sizeof(list*));

    // 更新坐标，把相关的点对加入树，然后把这个点加入网格
    for (uint32_t i = 0; i < point_num; i++) {
        refresh_hash(&object->point[i], N);
        distance_add(&point[i], hashmap, &distance, distance_list, N, point_num, cos_D);
        point_add(&point[i], hashmap, N);
    }

    const double slow_speed = 15.0;  // 这两个玄学，我也是瞎调的

    // 迭代的主循环，在这个循环以内的运算需要尽可能优化
    for (uint64_t i = 0; i < iteration; i++) {
        // 取出距离最近的一对点
        tree* nearest = (tree*)avl_tree_first(&distance);
        // 计算位移向量
        double move_rate = exp(i * (-slow_speed / iteration));
        vec3 move_vec_1 = nearest->point1->pos;
        sub(&move_vec_1, &nearest->point2->pos);
        scale(&move_vec_1, move_rate);
        vec3 move_vec_2 = move_vec_1;
        neg(&move_vec_2);
        // 移动这两个点，然后维护网格和树
        move_point(nearest->point1, &move_vec_1, hashmap, &distance, distance_list, N, point_num, cos_D);
        move_point(nearest->point2, &move_vec_2, hashmap, &distance, distance_list, N, point_num, cos_D);
    }

    object->angle = acos(((tree*)avl_tree_first(&distance))->cos_angle) * 180.0 / M_PI;

    free(hashmap);
    free(distance_list);
}

void creat_object(object* object, uint32_t seed, uint32_t point_num) {
    object->point_num = point_num;
    object->point = calloc(point_num, sizeof(list));
    dsfmt_init_gen_rand(&object->dsfmt, seed);
    for (uint32_t i = 0; i < point_num; i++) {
        sphere_point_picking(&object->point[i].pos, &object->dsfmt);
        object->point[i].id = i;
        object->point[i].next = NULL;
    }
}

#define OUTPUT_PRECISION "1.20"  //改这个可以修改输出的精度

void output_to_file(int3 version, time_t seed_0, object* object_list, uint32_t best_index, uint32_t point_num, uint64_t iteration, uint32_t repeat, double time) {
    char filename[256] = {0};
    sprintf(filename, "%u_%1.6lf.txt", point_num, object_list[best_index].angle);
    FILE* fp = fopen(filename, "w");

    if (fp == NULL)
        fp = stderr;

    fprintf(fp, "{");
    for (uint32_t i = 0; i < object_list[best_index].point_num; i++) {
        if (i < object_list[best_index].point_num - 1)
            fprintf(fp, "{%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf},", object_list[best_index].point[i].pos.x, object_list[best_index].point[i].pos.y, object_list[best_index].point[i].pos.z);
        else
            fprintf(fp, "{%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf,%" OUTPUT_PRECISION "lf}}\n", object_list[best_index].point[i].pos.x, object_list[best_index].point[i].pos.y, object_list[best_index].point[i].pos.z);
    }

    fprintf(fp, "version = %d.%d.%d\n", version.x, version.y, version.z);
    fprintf(fp, "seed_0 = %lu, point_num = %u, iteration = %lu, repeat = %u\n", seed_0, point_num, iteration, repeat);
    fprintf(fp, "best_index = %u, best_angle = %" OUTPUT_PRECISION "lf, time = %lf\n", best_index, object_list[best_index].angle, time);
    for (int i = 0; i < repeat; i++)
        fprintf(fp, "id = %4d , angle = %" OUTPUT_PRECISION "lf\n", i, object_list[i].angle);

    if (fp == stderr) {
        for (int i = 0; i < 3; i++)
            fprintf(stderr, "!!!写入文件时发生错误，尝试将结果输出至屏幕，请手动保存后再关闭程序!!!\n");
        while (1)
            getchar();
    }
    else {
        fclose(fp);
    }
}

int main(void) {
    const int3 version = {0, 2, 0};
    time_t seed_0 = time(NULL);

    uint32_t point_num;
    uint64_t iteration;
    uint32_t repeat;
    fprintf(stderr, "请依次输入节点数，优化迭代次数，重试次数，用空格分隔，然后按回车（示例：130 1000000 128）：\n");
    scanf("%" PRIu32 " %" PRIu64 " %" PRIu32, &point_num, &iteration, &repeat);

    double time_0 = omp_get_wtime();

    object* object_list = calloc(repeat, sizeof(object));

    int i;  // 如果写到for的括号里，msvc会报错
#pragma omp parallel for
    for (i = 0; i < repeat; i++) {
        uint32_t seed = (uint32_t)seed_0 + i;
        creat_object(&object_list[i], seed, point_num);
        tammes(&object_list[i], iteration);
        fprintf(stderr, "线程%4d优化完毕，最小夹角=%1.6lf\n", i, object_list[i].angle);
    }

    uint32_t best_index = 0;
    for (uint32_t i = 0; i < repeat; i++)
        if (object_list[i].angle > object_list[best_index].angle)
            best_index = i;

    double time_1 = omp_get_wtime();
    double time = time_1 - time_0;
    fprintf(stderr, "所有%" PRIu32 "个线程优化完毕，用时%lfs，被最大化的最小夹角=%1.16lf\n", repeat, time, object_list[best_index].angle);
    output_to_file(version, seed_0, object_list, best_index, point_num, iteration, repeat, time);

    for (uint32_t i = 0; i < repeat; i++)
        free(object_list[i].point);
    free(object_list);

    return 0;
}