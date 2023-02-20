#include "tammes.h"

#ifdef DEBUG
static int64_t adjust_count = 0;
#endif

typedef struct {
    int32_t point_idx1;
    int32_t point_idx2;
    double dot_of_point;
}point_pair_t;

int64_t equidistant_sum(int64_t an, int64_t am) {
    return ((an + am) * (am - an + 1)) >> 1;
}

void heap_push(point_pair_t point_pair, point_pair_t* point_pair_heap) {
    int64_t* count_ptr = (int64_t*)point_pair_heap;

    int64_t child_idx = ++(*count_ptr);
    int64_t father_idx = child_idx >> 1;

    while (child_idx > 2 && point_pair_heap[father_idx].dot_of_point < point_pair.dot_of_point) {
        point_pair_heap[child_idx] = point_pair_heap[father_idx];
        child_idx >>= 1;
        father_idx >>= 1;

#ifdef DEBUG
        adjust_count++;
#endif

    }

    point_pair_heap[child_idx] = point_pair;

}

double tammes(f64x3_t* point, int32_t point_num) {

    int64_t point_pair_num = equidistant_sum(1, point_num);
    point_pair_t* point_pair_heap = calloc(point_pair_num, sizeof(point_pair_t));

    for (int32_t i = 0;i < point_num;i++) {
        for (int32_t j = i + 1;j < point_num;j++) {
            double dot_of_point = f64x3_dot(&point[i], &point[j]);
            point_pair_t point_pair = { i, j, dot_of_point };
            heap_push(point_pair, point_pair_heap);
        }
    }

#ifdef DEBUG
    printf("heap count: %lld\n", *((int64_t*)point_pair_heap));
    printf("push count: %lld\n", adjust_count);
    printf("Max dot: %1.18lf\n", point_pair_heap[1].dot_of_point);
#endif

    return point_pair_heap[1].dot_of_point;

}