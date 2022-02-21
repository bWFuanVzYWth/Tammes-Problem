#pragma once

typedef struct {
    double x;
    double y;
    double z;
} vec3;

vec3 add(vec3 a, vec3 b);
vec3 sub(vec3 a, vec3 b);
vec3 mul(vec3 a, vec3 b);

vec3 scale(double a, vec3 b);
vec3 normalize(vec3 a);

double dot(vec3 a, vec3 b);

vec3 cross(vec3 a, vec3 b);