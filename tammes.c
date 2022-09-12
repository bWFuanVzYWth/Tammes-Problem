#include "tammes.h"

// 对于球面上的n点，总存在两个点，其距离<= D (Fejes Tóth, 1943)
double get_D(uint32_t point_num) {
    double tmp = 1.0 / sin((M_PI * point_num) / (6 * (point_num - 2)));
    return sqrt(4.0 - tmp * tmp);
}

inline uint64_t to_hash(vec3* pos, const uint32_t N) {
    int x = (int)((pos->x * 0.5 + 0.5) * N);
    int y = (int)((pos->y * 0.5 + 0.5) * N);
    int z = (int)((pos->z * 0.5 + 0.5) * N);
    return (x * N + y) * N + z;
}

// AVL树的比较函数，所有节点存放都在一个数组里，地址固定，因此地址的大小作为备选的排序方式
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

inline void point_add(list* point, list* hashmap) {
    point->last = &hashmap[point->hash];
    point->next = hashmap[point->hash].next;

    point->last->next = point;
    if (point->next != NULL)
        point->next->last = point;
}

inline void point_remove(list* point, list* hashmap) {
    point->last->next = point->next;
    if (point->next != NULL)
        point->next->last = point->last;
}

inline uint64_t to_p_tree(list* point1, list* point2, uint32_t point_num) {
    return point1->id * point_num + point2->id;
}

inline void distance_remove(list* point, list** hashmap_lut, struct avl_tree* distance, tree* distance_list, uint32_t point_num, const double cos_D) {
    uint64_t index = 28 * point->hash;
    while (hashmap_lut[index] != NULL) {
        list* p_del = hashmap_lut[index]->next;
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
        index++;
    }
}

inline void distance_add(list* point, list** hashmap_lut, struct avl_tree* distance, tree* distance_list, uint32_t point_num, const double cos_D) {
    uint64_t index = 28 * point->hash;
    while (hashmap_lut[index] != NULL) {
        list* p_add = hashmap_lut[index]->next;
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
        index++;
    }
}

inline void distance_refresh(list* point, list** hashmap_lut, struct avl_tree* distance, tree* distance_list, uint32_t point_num, const double cos_D, vec3* new_pos) {
    uint64_t index = 28 * point->hash;
    while (hashmap_lut[index] != NULL) {
        list* p_ref = hashmap_lut[index]->next;
        while (p_ref != NULL) {
            if (p_ref != point) {  // 跳过自己，不然会炸
                double old_cos_angle = dot(&point->pos, &p_ref->pos);
                double new_cos_angle = dot(new_pos, &p_ref->pos);
                int tmp_1 = isgreaterequal(old_cos_angle, cos_D);
                int tmp_2 = isgreaterequal(new_cos_angle, cos_D);

                if (tmp_1 || tmp_2) {
                    list* point1 = min(point, p_ref);
                    list* point2 = max(point, p_ref);
                    tree* p = distance_list + to_p_tree(point1, point2, point_num);
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
            }
            p_ref = p_ref->next;
        }
        index++;
    }
}

inline void move_point(list* point, vec3* new_pos, list* hashmap, list** hashmap_lut, struct avl_tree* distance, tree* distance_list, const uint32_t N, uint32_t point_num, const double cos_D) {
    uint64_t new_hash = to_hash(new_pos, N);
    if (point->hash == new_hash) {
        distance_refresh(point, hashmap_lut, distance, distance_list, point_num, cos_D, new_pos);
        point->pos = *new_pos;
    } else {
        point_remove(point, hashmap);
        distance_remove(point, hashmap_lut, distance, distance_list, point_num, cos_D);
        point->pos = *new_pos;
        point->hash = new_hash;
        distance_add(point, hashmap_lut, distance, distance_list, point_num, cos_D);
        point_add(point, hashmap);
    }
}

void tammes(object* object, uint64_t iteration) {
    list* point = object->point;
    list* best_point = object->best_point;
    uint32_t point_num = object->point_num;

    const double D = get_D(point_num);
    const uint32_t N = (uint32_t)floor(2.0 / D);
    const double L = 2.0 / N;
    const double cos_D = 1.0 - 0.5 * D * D - 1e-15;

    // 创建AVL树，然后创建距离列表，这个列表既是距离矩阵，也是AVL树的节点
    struct avl_tree distance;
    avl_tree_init(&distance, &tree_compare, sizeof(tree), AVL_OFFSET(tree, node));
    tree* distance_list = (tree*)calloc(point_num * point_num, sizeof(tree));
    // 创建空间划分网格，可以理解成hashmap，形式上是空的点的列表，仅使用指向下一个点的指针
    list* hashmap = (list*)calloc(N * N * N, sizeof(list));
    // 创建网格的邻接表（查找表），实际上是指向网格单元的指针数组，完成初始化后只读
    list** hashmap_lut = (list**)calloc(N * N * N * 28, sizeof(list*));

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
                            vec3 block_min = {block_min_x, block_min_y, block_min_z};
                            vec3 block_max = {block_max_x, block_max_y, block_max_z};
                            if (dot(&block_min, &block_min) < (1.0 + 1e-15) && dot(&block_max, &block_max) >= (1.0 - 1e-15))
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
        distance_add(point + i, hashmap_lut, &distance, distance_list, point_num, cos_D);
        point_add(point + i, hashmap);
    }

    // 一些玄学参数，自己看着调
    const double slow_down = 20.0;       // 控制移动的一个参数
    const double move_rate_0 = D / 2.0;  // 保证够大就行，似乎没啥区别

    object->angle = 2.0;

    // 迭代的主循环，在这个循环以内的运算需要尽可能优化
    for (uint64_t i = 0; i < iteration; i++) {
        // 取出距离最近的一对点
        tree* nearest = (tree*)avl_tree_first(&distance);

        if (nearest->cos_angle < object->angle) {
            memcpy(best_point, point, sizeof(list) * point_num);  // 记录最佳值
            object->angle = nearest->cos_angle;                   // 这里只是暂存一下余弦值
        }

        // 计算位移向量
        double move_rate = move_rate_0 * exp(i * (-slow_down / iteration));
        vec3 move_vec = nearest->point1->pos;
        sub(&move_vec, &nearest->point2->pos);
        double len = length(&move_vec);
        scale(&move_vec, move_rate / len);

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
    // 迭代的主循环结束

    object->angle = acos(object->angle) * 180.0 / M_PI;  // 这里才是把余弦值算成角度

    free(hashmap);
    free(distance_list);
}
