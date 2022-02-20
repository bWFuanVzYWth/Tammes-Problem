#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dSFMT.h"

#include "vec.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

vec3 sphere_point_picking(dsfmt_t* dsfmt) {
    // https://mathworld.wolfram.com/SpherePointPicking.html
    double x1, x2, tmp_1, tmp_2;
    do {
        x1 = dsfmt_genrand_open_open(dsfmt) * 2.0 - 1.0;  // (-1,1)
        x2 = dsfmt_genrand_open_open(dsfmt) * 2.0 - 1.0;  // (-1,1)
        tmp_1 = x1 * x1 + x2 * x2;
    } while (tmp_1 >= 1.0);
    tmp_2 = 2.0 * sqrt(1.0 - tmp_1);
    vec3 n = {x1 * tmp_2, x2 * tmp_2, 1.0 - 2.0 * tmp_1};
    return n;
}

typedef struct {
    vec3 point;
    size_t address;
} point_link;

typedef struct {
    int x;
    int y;
    int z;
} int3;

int3 hash_point(vec3 point, size_t n) {
    int3 hash;
    hash.x = (int)((point.x + 1.0) * 0.5 * n);
    hash.y = (int)((point.y + 1.0) * 0.5 * n);
    hash.z = (int)((point.z + 1.0) * 0.5 * n);
    return hash;
}

void find_nearest_point(vec3* point, size_t point_num, size_t* a, size_t* b) {
    // O(n)
    //对于球面上的n点，总存在两个点，其距离<=d (Fejes Tóth, 1943)
    //将空间以d为宽度划分为均匀网格，最小距离的两个点必然存在某3*3*3网格中
    double tmp_1 = 1.0 / sin((M_PI * point_num) / (6.0 * (point_num - 2)));
    double d = sqrt(4.0 - tmp_1 * tmp_1);
    size_t n = (size_t)floor(2.0 / d);

    size_t hash_map_size = sizeof(size_t) * n * n * n;
    size_t(*hash_map)[n][n] = malloc(hash_map_size);
    memset(hash_map, -1, hash_map_size);

    point_link* list = calloc(point_num, sizeof(point_link));

    double max_dot = -2.0;
    for (size_t i = 0; i < point_num; i++) {
        //计算需要搜索的网格的坐标范围
        int3 hash = hash_point(point[i], n);
        int3 hash_min = {max(hash.x - 1, 0), max(hash.y - 1, 0), max(hash.z - 1, 0)};
        int3 hash_max = {min(hash.x + 1, n - 1), min(hash.y + 1, n - 1), min(hash.z + 1, n - 1)};
        for (int x = hash_min.x; x <= hash_max.x; x++) {
            for (int y = hash_min.y; y <= hash_max.y; y++) {
                for (int z = hash_min.z; z <= hash_max.z; z++) {
                    // printf("aaaa\n");
                    if (hash_map[x][y][z] != -1) {
                        size_t address = hash_map[x][y][z];
                        do {
                            //计算距离，然后取出下一节链表的地址
                            double now_dot = dot(point[i], list[address].point);
                            // printf("now_dot=%lf\n", now_dot);
                            if (now_dot > max_dot) {
                                // max_dot = now_dot;
                                *a = address;
                                *b = i;
                            }
                            max_dot = fmax(max_dot, now_dot);

                            address = list[address].address;
                        } while (address != -1);
                    }
                }
            }
        }
        // printf("bbbb\n");

        if (hash_map[hash.x][hash.y][hash.z] == -1) {
            //加入链表并指向-1，加入hash map
            list[i].point = point[i];
            list[i].address = -1;
            hash_map[hash.x][hash.y][hash.z] = i;
        } else {
            //加入链表并指向原有的元素，修改hash map
            list[i].point = point[i];
            list[i].address = hash_map[hash.x][hash.y][hash.z];
            hash_map[hash.x][hash.y][hash.z] = i;
        }
    }
    // printf("max_dot=%1.20lf", max_dot);
    free(list);
    free(hash_map);
}

vec3 get_move_vec(double move_rate, dsfmt_t* dsfmt) {
    vec3 move_dir = sphere_point_picking(dsfmt);
    double move_len = move_rate * sqrt(log(1.0 / (1.0 - dsfmt_genrand_close_open(dsfmt))));
    return mul(move_len, move_dir);
}

void annealing(vec3* point, size_t point_num, dsfmt_t* dsfmt) {
    const size_t step = 1000000;
    const double cool_speed = 20.0;  //玄学，我也是瞎调的
    const double slow_speed = 10.0;  //玄学，我也是瞎调的

    double temperature = 1.0;
    double move_rate = 1.0;

    size_t a, b;
    find_nearest_point(point, point_num, &a, &b);
    double min_max_dot = dot(point[a], point[b]);

    for (size_t i = 0; i < step; i++) {
        temperature = exp(i * (-cool_speed / step));
        move_rate = exp(i * (-slow_speed / step));

        min_max_dot = dot(point[a], point[b]);
        vec3 new_point_a = normalize(add(point[a], get_move_vec(move_rate, dsfmt)));
        vec3 new_point_b = normalize(add(point[b], get_move_vec(move_rate, dsfmt)));

        double now_dot = dot(new_point_a, new_point_b);

        double dx = now_dot - min_max_dot;
        if (fmin(1.0, exp(dx / temperature)) < dsfmt_genrand_close_open(dsfmt)) {
            point[a] = new_point_a;
            point[b] = new_point_b;
            printf("step = %llu, min_angle = %lf\n", i, 180.0 * acos(min_max_dot) / M_PI);
        }
        find_nearest_point(point, point_num, &a, &b);
    }
}

vec3* tammes(size_t point_num, int32_t seed) {
    dsfmt_t dsfmt;
    dsfmt_init_gen_rand(&dsfmt, seed);

    vec3* point = (vec3*)calloc(point_num, sizeof(vec3));

    //随机生成点作为初始状态
    for (size_t i = 0; i < point_num; i++)
        point[i] = sphere_point_picking(&dsfmt);

    annealing(point, point_num, &dsfmt);

    return point;
}

void output_point(vec3* point, size_t point_num) {
    putchar('{');
    for (size_t i = 0; i < point_num; i++) {
        printf("{%1.6lf,%1.6lf,%1.6lf}", point[i].x, point[i].y, point[i].z);
        if (i < point_num - 1)
            putchar(',');
    }
    putchar('}');
}

int main(void) {
    int32_t seed = (int32_t)time(NULL);
    // int32_t seed = 1145141919;
    size_t point_num = 3000;
    vec3* point = tammes(point_num, seed);

    output_point(point, point_num);
    free(point);

    return 0;
}
