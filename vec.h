#pragma once

typedef struct {
    double x;
    double y;
    double z;
} vec3;

vec3 add(vec3 a, vec3 b);
vec3 mul(double a, vec3 b);

double dot(vec3 a, vec3 b);

vec3 normalize(vec3 a);