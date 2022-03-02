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
// size_t id;             // 4   点的编号
typedef struct list list;

struct list {
    list* next;
    vec3 pos;
    size_t hash;
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

void to_hash3(int3* hash3, vec3* pos, uint32_t N) {
    hash3->x = (int)((pos->x * 0.5 + 0.5) * N);
    hash3->y = (int)((pos->y * 0.5 + 0.5) * N);
    hash3->z = (int)((pos->z * 0.5 + 0.5) * N);
}

size_t to_hash(int x, int y, int z, const uint32_t N) {
    return (x * N + y) * N + z;
}

void refresh_hash_all(list* point, uint32_t N) {
    int3 hash3;
    to_hash3(&hash3, &point->pos, N);
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
    list* last = (list*)(hashmap + point->hash);
    list* here = hashmap[point->hash];
    while (point != here) {
        last = here;
        here = here->next;
    }
    last->next = here->next;
}

size_t to_p_tree(list* point1, list* point2, uint32_t point_num) {
    return point1->id * point_num + point2->id;
}

void distance_remove(list* point, list** hashmap, list*** hashmap_lut, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D) {
    size_t index = 27 * point->hash;
    for (int i = 0; i < 27 && hashmap_lut[index] != NULL; i++, index++) {
        list* p_del = *(hashmap_lut[index]);
        while (p_del != NULL) {
            double cos_angle = dot(&point->pos, &p_del->pos);
            if (isgreaterequal(cos_angle, cos_D)) {
                list* point1 = min(point, p_del);
                list* point2 = max(point, p_del);
                tree* p = distance_list + to_p_tree(point1, point2, point_num);
                avl_tree_remove(distance, p);
            }
            p_del = p_del->next;
        }
    }
}

void distance_add(list* point, list** hashmap, list*** hashmap_lut, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D) {
    size_t index = 27 * point->hash;
    for (int i = 0; i < 27 && hashmap_lut[index] != NULL; i++, index++) {
        list* p_add = *(hashmap_lut[index]);
        while (p_add != NULL) {
            double cos_angle = dot(&point->pos, &p_add->pos);
            if (isgreaterequal(cos_angle, cos_D)) {
                list* point1 = min(point, p_add);
                list* point2 = max(point, p_add);
                tree* p = distance_list + to_p_tree(point1, point2, point_num);
                p->point1 = point1;
                p->point2 = point2;
                p->cos_angle = cos_angle;
                avl_tree_add(distance, p);
            }
            p_add = p_add->next;
        }
    }
}

void distance_refresh(list* point, list** hashmap, list*** hashmap_lut, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D, vec3* new_pos) {
    size_t index = 27 * point->hash;
    for (int i = 0; i < 27 && hashmap_lut[index] != NULL; i++, index++) {
        list* p_ref = *(hashmap_lut[index]);
        while (p_ref != NULL) {
            if (p_ref != point) {  // 跳过自己
                double cos_angle = dot(&point->pos, &p_ref->pos);
                double new_cos_angle = dot(new_pos, &p_ref->pos);
                int tmp_1 = isgreaterequal(cos_angle, cos_D);
                int tmp_2 = isgreaterequal(new_cos_angle, cos_D);

                list* point1;
                list* point2;
                tree* p;
                if (tmp_1 || tmp_2) {
                    point1 = min(point, p_ref);
                    point2 = max(point, p_ref);
                    p = distance_list + to_p_tree(point1, point2, point_num);
                }
                if (tmp_1) {
                    avl_tree_remove(distance, p);
                }
                if (tmp_2) {
                    p->point1 = point1;
                    p->point2 = point2;
                    p->cos_angle = new_cos_angle;
                    avl_tree_add(distance, p);
                }
            }
            p_ref = p_ref->next;
        }
    }
}

void move_point(list* point, vec3* new_pos, list** hashmap, list*** hashmap_lut, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D) {
    int3 new_hash3;
    to_hash3(&new_hash3, new_pos, N);
    size_t new_hash = to_hash(new_hash3.x, new_hash3.y, new_hash3.z, N);

    if (point->hash == new_hash) {
        distance_refresh(point, hashmap, hashmap_lut, distance, distance_list, N, point_num, cos_D, new_pos);
        point->pos = *new_pos;
    } else {
        point_remove(point, hashmap, N);
        distance_remove(point, hashmap, hashmap_lut, distance, distance_list, N, point_num, cos_D);
        point->pos = *new_pos;
        point->hash = new_hash;
        distance_add(point, hashmap, hashmap_lut, distance, distance_list, N, point_num, cos_D);
        point_add(point, hashmap, N);
    }
}

void tammes(object* object, uint64_t iteration) {
    list* point = object->point;
    uint32_t point_num = object->point_num;

    const double D = get_D(point_num);
    const uint32_t N = (uint32_t)floor(2.0 / D);
    const double L = 2.0 / N;
    const double cos_D = 1.0 - 0.5 * D * D - 1e-7;

    // 初始化AVL树，然后创建距离矩阵，同时也是AVL树的节点
    struct avl_tree distance;
    avl_tree_init(&distance, &tree_compare, sizeof(tree), AVL_OFFSET(tree, node));
    tree* distance_list = (tree*)calloc(point_num * point_num, sizeof(tree));
    // 创建空间均匀划分网格，可以理解成hashmap
    list** hashmap = (list**)calloc(N * N * N, sizeof(list*));
    // 创建网格的邻域查找表
    list*** hashmap_lut = (list***)calloc(N * N * N * 27, sizeof(list**));

    // 对每一个网格，查找可能存在点的所有相邻网格，打表出奇迹
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                size_t index = 27 * to_hash(i, j, k, N);
                for (int x = max(i - 1, 0); x <= min(i + 1, N - 1); x++) {
                    for (int y = max(j - 1, 0); y <= min(j + 1, N - 1); y++) {
                        for (int z = max(k - 1, 0); z <= min(k + 1, N - 1); z++) {
                            int sign_x = x < N / 2 ? 0 : 1;
                            int sign_y = y < N / 2 ? 0 : 1;
                            int sign_z = z < N / 2 ? 0 : 1;
                            double block_max_x = (x + sign_x) * L - 1.0;
                            double block_max_y = (y + sign_y) * L - 1.0;
                            double block_max_z = (z + sign_z) * L - 1.0;
                            double block_min_x = (x + !sign_x) * L - 1.0;
                            double block_min_y = (y + !sign_y) * L - 1.0;
                            double block_min_z = (z + !sign_z) * L - 1.0;
                            vec3 block_min = {block_min_x, block_min_y, block_min_z};
                            vec3 block_max = {block_max_x, block_max_y, block_max_z};
                            if (dot(&block_min, &block_min) < 1.0 + 1e-7 && dot(&block_max, &block_max) >= 1.0 - 1e-7)
                                hashmap_lut[index++] = hashmap + to_hash(x, y, z, N);
                        }
                    }
                }
            }
        }
    }

    // 更新坐标，把相关的点对加入树，然后把这个点加入网格
    for (uint32_t i = 0; i < point_num; i++) {
        refresh_hash_all(object->point + i, N);
        distance_add(point + i, hashmap, hashmap_lut, &distance, distance_list, N, point_num, cos_D);
        point_add(point + i, hashmap, N);
    }

    const double slow_speed = 15.0;  // 这个玄学，我也是瞎调的

    // 迭代的主循环，在这个循环以内的运算需要尽可能优化
    for (uint64_t i = 0; i < iteration; i++) {
        // 取出距离最近的一对点
        tree* nearest = (tree*)avl_tree_first(&distance);
        // 计算位移向量
        double move_rate = exp(i * (-slow_speed / iteration));
        vec3 move_vec = nearest->point1->pos;
        sub(&move_vec, &nearest->point2->pos);
        scale(&move_vec, move_rate);
        //计算新的点坐标
        vec3 new_pos1 = nearest->point1->pos;
        add(&new_pos1, &move_vec);
        normalize(&new_pos1);
        vec3 new_pos2 = nearest->point2->pos;
        sub(&new_pos2, &move_vec);
        normalize(&new_pos2);

        // 移动这两个点，然后维护网格和树
        move_point(nearest->point1, &new_pos1, hashmap, hashmap_lut, &distance, distance_list, N, point_num, cos_D);
        move_point(nearest->point2, &new_pos2, hashmap, hashmap_lut, &distance, distance_list, N, point_num, cos_D);
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

#define OUTPUT_PRECISION "1.16"  //改这个可以修改输出的精度

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
    fprintf(fp, "seed_0 = %llu, point_num = %u, iteration = %llu, repeat = %u\n", seed_0, point_num, iteration, repeat);
    fprintf(fp, "best_index = %u, best_angle = %" OUTPUT_PRECISION "lf, time = %lf\n", best_index, object_list[best_index].angle, time);
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
        creat_object(object_list + i, seed, point_num);
        tammes(object_list + i, iteration);
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