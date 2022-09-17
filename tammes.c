#include "tammes.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

// struct point_t* next;     // 8   链表指针
// f64x3_t pos;              // 8*3 点的坐标
// uint64_t hash;           // 8   散列表
// uint32_t id;             // 4   点的编号
typedef struct point_t point_t;
typedef struct point_pair_t point_pair_t;

struct point_t {
    point_t* last;
    point_t* next;
    f64x3_t pos;
    uint64_t hash;
    uint32_t id;
};

// point_t* point1;          // 8     序号更小的点指针
// point_t* point2;          // 8     序号更大的点指针
// double cos_angle;      // 8     两点之间的夹角余弦
// struct avl_node node;  // 8*3+4 AVL树的指针
struct point_pair_t {
    point_t* point1;
    point_t* point2;
    double cos_angle;
    struct avl_node node;
};

// 对于球面上的n点，总存在两个点，其距离<= D (Fejes Tóth, 1943)
double get_D(uint32_t point_num) {
    double tmp = 1.0 / sin((M_PI * point_num) / (6 * (point_num - 2)));
    return sqrt(4.0 - tmp * tmp);
}

uint64_t to_hash(f64x3_t* pos, const uint32_t N) {
    int x = (int)((pos->x * 0.5 + 0.5) * N);
    int y = (int)((pos->y * 0.5 + 0.5) * N);
    int z = (int)((pos->z * 0.5 + 0.5) * N);
    return (x * N + y) * N + z;
}

// AVL树的比较函数，所有节点存放都在一个数组里，地址固定，因此地址的大小作为备选的排序方式
int tree_compare(const void* v_p1, const void* v_p2) {
    point_pair_t* p1 = (point_pair_t*)v_p1;
    point_pair_t* p2 = (point_pair_t*)v_p2;

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

void hashmap_add(point_t* point, point_t* hashmap) {
    point_t* hashmap_hash = hashmap + point->hash;
    point_t* hashmap_first = hashmap[point->hash].next;

    hashmap_hash->next = point;
    if (hashmap_first != NULL)
        hashmap_first->last = point;

    point->last = hashmap_hash;
    point->next = hashmap_first;
}

void hashmap_remove(point_t* point) {
    point->last->next = point->next;
    if (point->next != NULL)
        point->next->last = point->last;
}

point_pair_t* to_p_tree(point_pair_t* point_pair, point_t* point1, point_t* point2, const uint32_t point_num) {
    return point_pair + point1->id * point_num + point2->id;
}

void point_remove(point_t* point, point_t** hashmap_lut, struct avl_tree* distance, point_pair_t* point_pair, const uint32_t point_num, const double cos_D) {
    hashmap_remove(point);
    for (uint64_t index = 28 * point->hash; hashmap_lut[index] != NULL; index++) {
        for (point_t* p_del = hashmap_lut[index]->next; p_del != NULL; p_del = p_del->next) {
            double cos_angle = f64x3_dot(&point->pos, &p_del->pos);
            if (isless(cos_angle, cos_D))
                continue;
            point_t* point1 = min(point, p_del);
            point_t* point2 = max(point, p_del);
            point_pair_t* p = to_p_tree(point_pair, point1, point2, point_num);
            avl_tree_remove(distance, p);
        }
    }
}

void point_add(point_t* point, point_t* hashmap, point_t** hashmap_lut, struct avl_tree* distance, point_pair_t* point_pair, const uint32_t point_num, const double cos_D) {
    for (uint64_t index = 28 * point->hash; hashmap_lut[index] != NULL; index++) {
        for (point_t* p_add = hashmap_lut[index]->next; p_add != NULL; p_add = p_add->next) {
            double cos_angle = f64x3_dot(&point->pos, &p_add->pos);
            if (isless(cos_angle, cos_D))
                continue;
            point_t* point1 = min(point, p_add);
            point_t* point2 = max(point, p_add);
            point_pair_t* p = to_p_tree(point_pair, point1, point2, point_num);
            p->point1 = point1;
            p->point2 = point2;
            p->cos_angle = cos_angle;
            avl_tree_add(distance, p);
        }
    }
    hashmap_add(point, hashmap);
}

void point_refresh(point_t* point, point_t** hashmap_lut, struct avl_tree* distance, point_pair_t* point_pair, const uint32_t point_num, const double cos_D, f64x3_t* new_pos) {
    for (uint64_t index = 28 * point->hash; hashmap_lut[index] != NULL; index++) {
        for (point_t* p_ref = hashmap_lut[index]->next; p_ref != NULL; p_ref = p_ref->next) {
            if (p_ref != point) {  // 跳过自己，不然会炸
                double old_cos_angle = f64x3_dot(&point->pos, &p_ref->pos);
                double new_cos_angle = f64x3_dot(new_pos, &p_ref->pos);
                int tmp_1 = isless(old_cos_angle, cos_D);
                int tmp_2 = isless(new_cos_angle, cos_D);

                if (tmp_1 && tmp_2)
                    continue;
                point_t* point1 = min(point, p_ref);
                point_t* point2 = max(point, p_ref);
                point_pair_t* p = to_p_tree(point_pair, point1, point2, point_num);
                if (!tmp_1) {
                    avl_tree_remove(distance, p);
                }
                if (!tmp_2) {
                    p->point1 = point1;
                    p->point2 = point2;
                    p->cos_angle = new_cos_angle;
                    avl_tree_add(distance, p);
                }
            }
        }
    }
}

void move_point(point_t* point, f64x3_t* new_pos, point_t* hashmap, point_t** hashmap_lut, struct avl_tree* distance, point_pair_t* point_pair, const uint32_t N, const uint32_t point_num, const double cos_D) {
    uint64_t new_hash = to_hash(new_pos, N);
    if (point->hash == new_hash) {
        point_refresh(point, hashmap_lut, distance, point_pair, point_num, cos_D, new_pos);
        point->pos = *new_pos;
    } else {
        point_remove(point, hashmap_lut, distance, point_pair, point_num, cos_D);
        point->pos = *new_pos;
        point->hash = new_hash;
        point_add(point, hashmap, hashmap_lut, distance, point_pair, point_num, cos_D);
    }
}

void pos_to_point(f64x3_t* pos, point_t* point, const uint32_t point_num) {
    for (int i = 0; i < point_num; i++) {
        point[i].id = i;
        point[i].pos.x = pos[i].x;
        point[i].pos.y = pos[i].y;
        point[i].pos.z = pos[i].z;
    }
}

void point_to_pos(point_t* point, f64x3_t* pos, const uint32_t point_num) {
    for (int i = 0; i < point_num; i++) {
        pos[i].x = point[i].pos.x;
        pos[i].y = point[i].pos.y;
        pos[i].z = point[i].pos.z;
    }
}

double tammes(f64x3_t* pos, const uint32_t point_num, const uint64_t iteration) {
    point_t* point = calloc(point_num, sizeof(point_t));
    pos_to_point(pos, point, point_num);

    const double D = get_D(point_num);
    const uint32_t N = (uint32_t)floor(2.0 / D);
    const double L = 2.0 / N;
    const double cos_D = 1.0 - 0.5 * D * D;

    // 创建AVL树，然后创建距离列表，这个列表既是距离矩阵，也是AVL树的节点
    struct avl_tree distance;
    avl_tree_init(&distance, &tree_compare, sizeof(point_pair_t), AVL_OFFSET(point_pair_t, node));
    point_pair_t* point_pair = (point_pair_t*)calloc(point_num * point_num, sizeof(point_pair_t));
    // 创建空间划分网格，可以理解成hashmap，形式上是空的点的列表，仅使用指向下一个点的指针
    point_t* hashmap = (point_t*)calloc(N * N * N, sizeof(point_t));
    // 创建网格的邻接表（查找表），实际上是指向网格单元的指针数组，完成初始化后只读
    point_t** hashmap_lut = (point_t**)calloc(N * N * N * 28, sizeof(point_t*));

    // 初始化网格的邻接表，每个网格最多与3*3*3=27个其他网格相接
    // 遍历每一个网格
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                uint64_t index = 28 * ((i * N + j) * N + k);
                //查找周围3*3*3网格
                for (int x = max(i - 1, 0); x <= min(i + 1, N - 1); x++) {
                    for (int y = max(j - 1, 0); y <= min(j + 1, N - 1); y++) {
                        for (int z = max(k - 1, 0); z <= min(k + 1, N - 1); z++) {
                            //检查这个网格是否与球面相交，相交说明有可能存在点，打表
                            int sign_x = x < N / 2 ? 0 : 1;
                            int sign_y = y < N / 2 ? 0 : 1;
                            int sign_z = z < N / 2 ? 0 : 1;
                            double block_max_x = (x + sign_x) * L - 1.0;
                            double block_max_y = (y + sign_y) * L - 1.0;
                            double block_max_z = (z + sign_z) * L - 1.0;
                            double block_min_x = (x + !sign_x) * L - 1.0;
                            double block_min_y = (y + !sign_y) * L - 1.0;
                            double block_min_z = (z + !sign_z) * L - 1.0;
                            f64x3_t block_min = {block_min_x, block_min_y, block_min_z};
                            f64x3_t block_max = {block_max_x, block_max_y, block_max_z};
                            if (f64x3_dot(&block_min, &block_min) < 1.0 && f64x3_dot(&block_max, &block_max) >= 1.0)
                                hashmap_lut[index++] = &hashmap[(x * N + y) * N + z];
                        }
                    }
                }
            }
        }
    }

    // 初始化AVL树，遍历所有节点，记录这个节点和所有相近的节点的距离，并根据坐标放进空间划分网格
    for (uint32_t i = 0; i < point_num; i++) {
        point[i].hash = to_hash(&point[i].pos, N);
        point_add(point + i, hashmap, hashmap_lut, &distance, point_pair, point_num, cos_D);
    }

    // 一些玄学参数，自己看着调
    const double slow_down = 20.0;       // 控制移动的一个参数
    const double move_rate_0 = D / 2.0;  // 保证够大就行，似乎没啥区别

    double max_cos = 1.0;

    // 迭代的主循环，在这个循环以内的运算需要尽可能优化
    for (uint64_t i = 0; i < iteration; i++) {
        // 取出距离最近的一对点
        point_pair_t* nearest = (point_pair_t*)avl_tree_first(&distance);

        if (nearest->cos_angle < max_cos) {
            point_to_pos(point, pos, point_num);
            max_cos = nearest->cos_angle;
        }

        // 计算位移向量
        double move_rate = move_rate_0 * exp(i * (-slow_down / iteration));
        f64x3_t move_vec = nearest->point1->pos;
        f64x3_sub(&move_vec, &nearest->point2->pos);
        double len = f64x3_length(&move_vec);
        f64x3_mul(&move_vec, move_rate / len);

        //计算新的点坐标
        f64x3_t new_pos1 = nearest->point1->pos;
        f64x3_add(&new_pos1, &move_vec);
        f64x3_normalize(&new_pos1);
        f64x3_t new_pos2 = nearest->point2->pos;
        f64x3_sub(&new_pos2, &move_vec);
        f64x3_normalize(&new_pos2);

        // 移动这两个点，然后维护网格和树
        move_point(nearest->point1, &new_pos1, hashmap, hashmap_lut, &distance, point_pair, N, point_num, cos_D);
        move_point(nearest->point2, &new_pos2, hashmap, hashmap_lut, &distance, point_pair, N, point_num, cos_D);
    }
    // 迭代的主循环结束

    double angle = acos(max_cos) * 180.0 / M_PI;

    free(hashmap);
    free(point_pair);

    return angle;
}