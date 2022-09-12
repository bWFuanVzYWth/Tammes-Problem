#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "vec.h"
#include "avlmini.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

// list 是节点的坐标列表的最小单位，以双向链表的形式储存
// struct list* next;     // 8   链表指针
// vec3 pos;              // 8*3 点的坐标
// uint64_t hash;           // 8   散列表
// uint32_t id;             // 4   点的编号
typedef struct list list;

struct list {
    list* last;
    list* next;
    vec3 pos;
    uint64_t hash;
    uint32_t id;
};

// object 是一个待优化的对象，可以创建多个object同时优化
// list* point;         // 32*n   点的坐标列表
// double angle;        // 8      最小距离点对的夹角
// uint32_t point_num;  // 4      点的数量
typedef struct {
    list* point;
    list* best_point;
    double angle;
    uint32_t point_num;
} object;

// tree 是优化算法中用到的数据结构，维护了不同节点之间的距离
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



void tammes(object* object, uint64_t iteration);
