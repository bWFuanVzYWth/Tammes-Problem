#include "tammes.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

typedef struct point_t point_t;
typedef struct point_pair_t point_pair_t;

// f64x3_t pos;             // 8*3 点的坐标
// struct point_t* last;    // 8   链表指针
// struct point_t* next;    // 8   链表指针
// uint64_t hash;           // 8   散列表
// uint32_t id;             // 4   点的编号
struct point_t {
    f64x3_t pos;
    point_t* last;
    point_t* next;
    uint64_t hash;
    uint32_t id;
};

// struct avl_node node;    // 8*3+4 AVL树的指针
// point_t* point1;         // 8     序号更小的点指针
// point_t* point2;         // 8     序号更大的点指针
// double cos;              // 8     两点之间的夹角余弦
struct point_pair_t {
    struct avl_node node;
    point_t* point1;
    point_t* point2;
    double cos;
};

double cos_to_deg(double cos) {
    return acos(cos) * (180.0 / M_PI);
}

double cos_to_dis(double cos) {
    return sqrt(2.0 - 2.0 * cos);
}

double dis_to_cos(double dis) {
    return 1.0 - 0.5 * dis * dis;
}

// AVL树的比较函数，优先比较点对间的距离，距离相等时比较指针的大小（极小概率）
int tree_compare(const void* v_p1, const void* v_p2) {
    point_pair_t* p1 = (point_pair_t*)v_p1;
    point_pair_t* p2 = (point_pair_t*)v_p2;

    // 先检查是否为同一节点，避免浮点误差
    if (p1 == p2)
        return 0;
    // 比较点对间的距离，cos越大夹角越小所以是反的
    if (isless(p2->cos, p1->cos))
        return -1;
    if (isless(p1->cos, p2->cos))
        return 1;
    // 距离相等时比较指针的大小（极小概率）
    if (p1 < p2)
        return -1;
    else
        return 1;
}

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

void point_remove(point_t* point, point_t** hashmap_lut, struct avl_tree* point_pair_tree, point_pair_t* point_pair, const uint32_t point_num, const double cos_D) {
    hashmap_remove(point);
    for (uint64_t index = 28 * point->hash; hashmap_lut[index] != NULL; index++) {
        for (point_t* p_del = hashmap_lut[index]->next; p_del != NULL; p_del = p_del->next) {
            double cos = f64x3_dot(&point->pos, &p_del->pos);
            if (isless(cos, cos_D))
                continue;
            point_t* point1 = min(point, p_del);
            point_t* point2 = max(point, p_del);
            point_pair_t* p = to_p_tree(point_pair, point1, point2, point_num);
            avl_tree_remove(point_pair_tree, p);
        }
    }
}

void point_add(point_t* point, point_t* hashmap, point_t** hashmap_lut, struct avl_tree* point_pair_tree, point_pair_t* point_pair, const uint32_t point_num, const double cos_D) {
    for (uint64_t index = 28 * point->hash; hashmap_lut[index] != NULL; index++) {
        for (point_t* p_add = hashmap_lut[index]->next; p_add != NULL; p_add = p_add->next) {
            double cos = f64x3_dot(&point->pos, &p_add->pos);
            if (isless(cos, cos_D))
                continue;
            point_t* point1 = min(point, p_add);
            point_t* point2 = max(point, p_add);
            point_pair_t* p = to_p_tree(point_pair, point1, point2, point_num);
            p->point1 = point1;
            p->point2 = point2;
            p->cos = cos;
            avl_tree_add(point_pair_tree, p);
        }
    }
    hashmap_add(point, hashmap);
}

void point_refresh(point_t* point, point_t** hashmap_lut, struct avl_tree* point_pair_tree, point_pair_t* point_pair, const uint32_t point_num, const double cos_D, f64x3_t* new_pos) {
    for (uint64_t index = 28 * point->hash; hashmap_lut[index] != NULL; index++) {
        for (point_t* p_ref = hashmap_lut[index]->next; p_ref != NULL; p_ref = p_ref->next) {
            if (p_ref == point)
                continue;
            double old_cos = f64x3_dot(&point->pos, &p_ref->pos);
            double new_cos = f64x3_dot(new_pos, &p_ref->pos);
            int tmp_1 = isless(old_cos, cos_D);
            int tmp_2 = isless(new_cos, cos_D);

            if (tmp_1 && tmp_2)
                continue;
            point_t* point1 = min(point, p_ref);
            point_t* point2 = max(point, p_ref);
            point_pair_t* p = to_p_tree(point_pair, point1, point2, point_num);
            if (!tmp_1) {
                avl_tree_remove(point_pair_tree, p);
            }
            if (!tmp_2) {
                p->point1 = point1;
                p->point2 = point2;
                p->cos = new_cos;
                avl_tree_add(point_pair_tree, p);
            }
        }
    }
}

void move_point(point_t* point, f64x3_t* new_pos, point_t* hashmap, point_t** hashmap_lut, struct avl_tree* point_pair_tree, point_pair_t* point_pair, const uint32_t N, const uint32_t point_num, const double cos_D) {
    uint64_t new_hash = to_hash(new_pos, N);
    if (point->hash == new_hash) {
        point_refresh(point, hashmap_lut, point_pair_tree, point_pair, point_num, cos_D, new_pos);
        point->pos = *new_pos;
    } else {
        point_remove(point, hashmap_lut, point_pair_tree, point_pair, point_num, cos_D);
        point->pos = *new_pos;
        point->hash = new_hash;
        point_add(point, hashmap, hashmap_lut, point_pair_tree, point_pair, point_num, cos_D);
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

void check_333(point_t** hashmap_lut, point_t* hashmap, int i, int j, int k, const uint32_t N, const double L) {
    uint64_t index = 28 * ((i * N + j) * N + k);
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
                f64x3_t block_min = {block_min_x, block_min_y, block_min_z};
                f64x3_t block_max = {block_max_x, block_max_y, block_max_z};
                if (f64x3_dot(&block_min, &block_min) < 1.0 && f64x3_dot(&block_max, &block_max) >= 1.0)
                    hashmap_lut[index++] = &hashmap[(x * N + y) * N + z];
            }
        }
    }
}

void hashmap_init(point_t** hashmap_lut, point_t* hashmap, const uint32_t N, const double L) {
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < N; k++)
                check_333(hashmap_lut, hashmap, i, j, k, N, L);
}

double tammes(f64x3_t* pos, const uint32_t point_num, const uint64_t iteration, const uint32_t mode) {
    const double D = get_D(point_num);
    const double cos_D = dis_to_cos(D);
    const uint32_t N = (uint32_t)floor(2.0 / D);
    const double L = 2.0 / N;

    point_t* point = calloc(point_num, sizeof(point_t));
    point_pair_t* point_pair = (point_pair_t*)calloc(point_num * point_num, sizeof(point_pair_t));
    point_t* hashmap = (point_t*)calloc(N * N * N, sizeof(point_t));
    point_t** hashmap_lut = (point_t**)calloc(N * N * N * 28, sizeof(point_t*));

    pos_to_point(pos, point, point_num);
    hashmap_init(hashmap_lut, hashmap, N, L);
    struct avl_tree point_pair_tree;
    avl_tree_init(&point_pair_tree, &tree_compare, sizeof(point_pair_t), AVL_OFFSET(point_pair_t, node));
    for (uint32_t i = 0; i < point_num; i++) {
        point[i].hash = to_hash(&point[i].pos, N);
        point_add(point + i, hashmap, hashmap_lut, &point_pair_tree, point_pair, point_num, cos_D);
    }

    // 一些玄学参数，自己看着调
    const double slow_down = 20.0;
    const double move_rate_0 = D / 2.0;

    double max_cos = 1.0;

    // 迭代的主循环，在这个循环以内尽量优化
    for (uint64_t i = 0; i < iteration; i++) {
        point_pair_t* nearest = (point_pair_t*)avl_tree_first(&point_pair_tree);

        if (nearest->cos < max_cos) {
            max_cos = nearest->cos;
            point_to_pos(point, pos, point_num);
        }

        double move_rate = move_rate_0 * exp(i * (-slow_down / iteration));

        f64x3_t move_vec = nearest->point1->pos;
        f64x3_sub(&move_vec, &nearest->point2->pos);
        double len = f64x3_length(&move_vec);
        f64x3_mul(&move_vec, move_rate / len);

        f64x3_t new_pos1 = nearest->point1->pos;
        f64x3_add(&new_pos1, &move_vec);
        f64x3_normalize(&new_pos1);
        f64x3_t new_pos2 = nearest->point2->pos;
        f64x3_sub(&new_pos2, &move_vec);
        f64x3_normalize(&new_pos2);

        move_point(nearest->point1, &new_pos1, hashmap, hashmap_lut, &point_pair_tree, point_pair, N, point_num, cos_D);
        move_point(nearest->point2, &new_pos2, hashmap, hashmap_lut, &point_pair_tree, point_pair, N, point_num, cos_D);
    }
    // 迭代的主循环结束

    free(point);
    free(point_pair);
    free(hashmap);
    free(hashmap_lut);

    return cos_to_deg(max_cos);
}