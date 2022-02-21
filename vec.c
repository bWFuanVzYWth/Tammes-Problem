#include <math.h>

#include "vec.h"

vec3 add(vec3 a, vec3 b) {
    vec3 c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    c.z = a.z + b.z;
    return c;
}

vec3 sub(vec3 a, vec3 b) {
    vec3 c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    c.z = a.z - b.z;
    return c;
}

vec3 mul(vec3 a, vec3 b) {
    vec3 c;
    c.x = a.x * b.x;
    c.y = a.y * b.y;
    c.z = a.z * b.z;
    return c;
}

vec3 scale(double a, vec3 b) {
    vec3 c;
    c.x = a * b.x;
    c.y = a * b.y;
    c.z = a * b.z;
    return c;
}

vec3 normalize(vec3 a) {
    double rcp_length = 1.0 / sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
    vec3 b;
    b.x = a.x * rcp_length;
    b.y = a.y * rcp_length;
    b.z = a.z * rcp_length;
    return b;
}

double dot(vec3 a, vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec3 cross(vec3 a, vec3 b) {
    vec3 c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}
