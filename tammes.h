#ifndef _TAMMES_HEAD_
#define _TAMMES_HEAD_

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "avlmini.h"
#include "vec.h"

double tammes(f64x3_t* pos, uint32_t point_num, uint64_t iteration, uint32_t mode);

#endif