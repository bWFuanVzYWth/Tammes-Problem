#define _USE_MATH_DEFINES
#include <inttypes.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "avlmini.h"
#include "dSFMT.h"
#include "vec.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

// struct list* next;  // 8   链表指针
// vec3 pos;           // 8*3 点的坐标
// int3 hash3;         // 4*3 空间划分网格
// int hash;           // 4   散列表
// size_t id;          // 4   点的编号
struct list {
    struct list* next;
    vec3 pos;
    int3 hash3;
    int hash;
    uint32_t id;
};
typedef struct list list;

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
// struct avl_node node;  // 8*3+4 AVL树
typedef struct {
    list* point1;
    list* point2;
    double cos_angle;
    struct avl_node node;
} tree;

// 输出坐标列表
void output_point(object* object) {
    putchar('\n');
    putchar('{');
    for (uint32_t i = 0; i < object->point_num; i++) {
        printf("{%1.6lf,%1.6lf,%1.6lf}", object->point[i].pos.x, object->point[i].pos.y, object->point[i].pos.z);
        if (i < object->point_num - 1)
            putchar(',');
    }
    putchar('}');
    putchar('\n');
}

// 在单位圆内随机生成一个点
void circle_point_picking(double* x1, double* x2, double* r2, dsfmt_t* dsfmt) {
    do {
        *x1 = dsfmt_genrand_open_open(dsfmt) * 2.0 - 1.0;
        *x2 = dsfmt_genrand_open_open(dsfmt) * 2.0 - 1.0;
        *r2 = *x1 * *x1 + *x2 * *x2;
    } while (*r2 >= 1.0);
}

// 在单位球面上随机生成一个点
void sphere_point_picking(vec3* v, dsfmt_t* dsfmt) {
    double x1, x2, r2, tmp;
    circle_point_picking(&x1, &x2, &r2, dsfmt);
    tmp = 2.0 * sqrt(1.0 - r2);
    v->x = x1 * tmp;
    v->y = x2 * tmp;
    v->z = 1.0 - 2.0 * r2;
}

// AVL树的比较函数
// 所有节点存放都在一个数组里，每个节点的地址固定，因此地址的大小也作为比较方式之一
int tree_compare(const void* v_p1, const void* v_p2) {
    tree* p1 = (tree*)v_p1;
    tree* p2 = (tree*)v_p2;

    // 先检查是否为同一节点，避免浮点误差
    if (p1->point1 == p2->point1 && p1->point2 == p2->point2)
        return 0;

    // 距离不同时的比较，使用math.h中的函数比较浮点数大小，cos越大夹角越小所以是反的
    if (isgreater(p1->cos_angle, p2->cos_angle))
        return -1;
    if (isless(p1->cos_angle, p2->cos_angle))
        return 1;

    // 相等时的处理，应该很少走到这里
    if ((p1->point1 < p2->point1) || (p1->point1 == p2->point1 && p1->point2 < p2->point2))
        return -1;
    else
        return 1;
}

// 将空间坐标映射到空间划分网格坐标
void to_hash3(int3* hash3, vec3* point, uint32_t N) {
    hash3->x = (int)((point->x * 0.5 + 0.5) * N);
    hash3->y = (int)((point->y * 0.5 + 0.5) * N);
    hash3->z = (int)((point->z * 0.5 + 0.5) * N);
}

// 将空间划分网格坐标映射到hash
int to_hash(int3* hash3, uint32_t N) {
    return (hash3->x * N + hash3->y) * N + hash3->z;
}

void refresh_hash(list* point, uint32_t N) {
    to_hash3(&point->hash3, &point->pos, N);
    point->hash = to_hash(&point->hash3, N);
}

void hashmap_add(list* point, list** hashmap, const uint32_t N) {
    point->next = hashmap[point->hash];
    hashmap[point->hash] = point;
}

void hashmap_remove(list* point, list** hashmap, const uint32_t N) {
    list* last = (list*)&hashmap[point->hash];
    list* here = hashmap[point->hash];

    // 向后查找链表找到需要删除的节点
    while (point != here) {
        last = here;
        here = (list*)here->next;
    }

    last->next = here->next;
}

// 对于球面上的n点，总存在两个点，其距离<=D (Fejes Tóth, 1943)
double get_D(uint32_t point_num) {
    double tmp = 1.0 / sin((M_PI * point_num) / (6.0 * (point_num - 2)));
    return sqrt(4.0 - tmp * tmp);
}

void move_point(list* point, vec3* move_vec, list** hashmap, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D) {
    // 移除hashmap
    hashmap_remove(point, hashmap, N);

    // 移除旧坐标附近的所有连接
    int3 hash3_del_min = {max(point->hash3.x - 1, 0), max(point->hash3.y - 1, 0), max(point->hash3.z - 1, 0)};
    int3 hash3_del_max = {min(point->hash3.x + 1, N - 1), min(point->hash3.y + 1, N - 1), min(point->hash3.z + 1, N - 1)};
    int3 hash3_del;
    for (hash3_del.x = hash3_del_min.x; hash3_del.x <= hash3_del_max.x; hash3_del.x++) {
        for (hash3_del.y = hash3_del_min.y; hash3_del.y <= hash3_del_max.y; hash3_del.y++) {
            for (hash3_del.z = hash3_del_min.z; hash3_del.z <= hash3_del_max.z; hash3_del.z++) {
                size_t del_index = to_hash(&hash3_del, N);
                list* p_del = hashmap[del_index];
                while (p_del != NULL) {
                    double cos_angle = dot(&point->pos, &p_del->pos);
                    if (isgreaterequal(cos_angle, cos_D)) {
                        tree del;
                        del.point1 = min(point, p_del);
                        del.point2 = max(point, p_del);
                        del.cos_angle = cos_angle;
                        tree* p_del = avl_tree_find(distance, &del);
                        avl_tree_remove(distance, p_del);
                    }
                    p_del = (list*)p_del->next;
                }
            }
        }
    }

    // 移动点
    add(&point->pos, move_vec);
    normalize(&point->pos);
    refresh_hash(point, N);

    // 加入新坐标附近的所有连接
    int3 hash3_add_min = {max(point->hash3.x - 1, 0), max(point->hash3.y - 1, 0), max(point->hash3.z - 1, 0)};
    int3 hash3_add_max = {min(point->hash3.x + 1, N - 1), min(point->hash3.y + 1, N - 1), min(point->hash3.z + 1, N - 1)};
    int3 hash3_add;
    for (hash3_add.x = hash3_add_min.x; hash3_add.x <= hash3_add_max.x; hash3_add.x++) {
        for (hash3_add.y = hash3_add_min.y; hash3_add.y <= hash3_add_max.y; hash3_add.y++) {
            for (hash3_add.z = hash3_add_min.z; hash3_add.z <= hash3_add_max.z; hash3_add.z++) {
                size_t add_index = to_hash(&hash3_add, N);
                list* p_add = hashmap[add_index];
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
                    p_add = (list*)p_add->next;
                }
            }
        }
    }

    // 加入hashmap
    hashmap_add(point, hashmap, N);
}

// 传入一个对象，通过算法优化
void tammes(object* object, uint64_t iteration) {
    list* point = object->point;
    uint32_t point_num = object->point_num;

    const double D = get_D(point_num);  //对于球面上的n点，总存在两个点，其距离<=d (Fejes Tóth, 1943)
    const uint32_t N = (uint32_t)floor(2.0 / D);
    const double cos_D = 1.0 - 0.5 * D * D - 1e-15;  // 余弦定理，-1e-15保证浮点精度

    // 创建空间均匀划分网格
    list** hashmap = (list**)calloc(N * N * N, sizeof(list*));
    for (uint32_t i = 0; i < point_num; i++) {
        refresh_hash(&object->point[i], N);
        hashmap_add(&point[i], hashmap, N);
    }

    // 创建AVL树的树根
    struct avl_tree distance;
    avl_tree_init(&distance, &tree_compare, sizeof(tree), AVL_OFFSET(tree, node));
    // 创建AVL树的节点
    const size_t distance_num = point_num * point_num;
    tree* distance_list = (tree*)calloc(distance_num, sizeof(tree));

    // 遍历两个点的组合，把所有点的连接关系存进树里，其实这里也能用空间划分网格但是懒，<1ms
    for (uint32_t i = 0; i < point_num - 1; i++) {
        for (uint32_t j = i + 1; j < point_num; j++) {
            double cos_angle = dot(&point[i].pos, &point[j].pos);
            if (isgreaterequal(cos_angle, cos_D)) {
                tree* p = &distance_list[i * point_num + j];
                p->point1 = &point[i];
                p->point2 = &point[j];
                p->cos_angle = cos_angle;
                avl_tree_add(&distance, p);
            }
        }
    }

    const double slow_speed = 15.0;  // 这两个玄学，我也是瞎调的
    double move_rate = 1.0;

    // 迭代的主循环 TODO 哪怕一点点优化
    for (uint64_t i = 0; i < iteration; i++) {
        // 取出距离最近的一对点
        tree* nearest = (tree*)avl_tree_first(&distance);
        vec3 pos1 = nearest->point1->pos;
        vec3 pos2 = nearest->point2->pos;

        // 计算位移向量（已经变成汇编的形状了
        move_rate = exp(i * (-slow_speed / iteration));
        vec3 move_vec_1 = pos1;
        sub(&move_vec_1, &pos2);
        scale(&move_vec_1, move_rate);
        vec3 move_vec_2 = move_vec_1;
        neg(&move_vec_2);

        move_point(nearest->point1, &move_vec_1, hashmap, &distance, distance_list, N, point_num, cos_D);
        move_point(nearest->point2, &move_vec_2, hashmap, &distance, distance_list, N, point_num, cos_D);
    }
    object->angle = acos(((tree*)avl_tree_first(&distance))->cos_angle) * 180.0 / M_PI;

    free(hashmap);
    free(distance_list);
}

// 写入参数、申请内存、设置伪随机种子、生成点的初始排列
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

int main(void) {
    uint32_t seed_0 = (uint32_t)time(NULL);

    uint32_t point_num;
    uint64_t iteration;
    uint32_t repeat;
    fprintf(stderr, "提示：可在打开软件前通过>>符号将输出重定向至文件而不是屏幕\n");
    fprintf(stderr, "请依次输入节点数，优化迭代次数，重试次数，用空格分隔，然后按回车（示例：130 1000000 128）：\n");
    scanf("%" PRIu32 " %" PRIu64 " %" PRIu32, &point_num, &iteration, &repeat);

    object* object_list = calloc(repeat, sizeof(object));

#pragma omp parallel for
    for (uint32_t i = 0; i < repeat; i++) {
        uint32_t seed = seed_0 + i;
        creat_object(&object_list[i], seed, point_num);
        tammes(&object_list[i], iteration);
    }

    uint32_t best_index = 0;
    for (uint32_t i = 0; i < repeat; i++)
        if (object_list[i].angle > object_list[best_index].angle)
            best_index = i;

    fprintf(stderr, "所有%" PRIu32 "个线程优化完毕，被最大化的最小夹角 = %lf\n", repeat, object_list[best_index].angle);
    output_point(&object_list[best_index]);

    for (uint32_t i = 0; i < repeat; i++)
        free(object_list[i].point);
    free(object_list);

    getchar();

    return 0;
}