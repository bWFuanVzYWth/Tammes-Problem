#ifndef _TAMMES
#define _TAMMES

#include <stdlib.h>

#ifdef DEBUG
#include <stdio.h>
#endif

#include "vec.h"

// 优化输入的点的坐标，并返回优化后点的最小夹角
double tammes(f64x3_t *point, int32_t point_num);

#endif