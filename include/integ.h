#ifndef _INTEG_H
#define _INTEG_H

#include <math.h>
#include <stdbool.h>

#include "matrix.h"
#include "vector.h"

double integ_trap(vector_t *y, vector_t *x);

typedef struct {
    double (*f)(double);
    double x0;
    double xf;
    unsigned n;
    matrix_t *R;
} romberg_t;

romberg_t *new_romberg(double (*f)(double), double x0, double xf, unsigned n);

void free_romberg(romberg_t *romb);

void print_romberg(romberg_t *romb);

void solver_romberg(romberg_t *romb);

#endif
