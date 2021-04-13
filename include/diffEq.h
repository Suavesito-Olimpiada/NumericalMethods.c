#ifndef DIFFEQ_H
#define DIFFEQ_H 1

#include <stdlib.h>

#include "vector.h"

/* Initial value problem */
typedef struct {
    double y0;
    double (*f)(double, double);
} pvi_t;

/* create initial value problem */
pvi_t *new_pvi(double (*f)(double, double), double y0);

/* delete initial value problem */
void free_pvi(pvi_t *pvi);

/* print initil value problem */
void print_pvi(pvi_t *pvi);

/* Runge-Kutta 2nd order problem */
typedef struct {
    double a;
    double b;
    unsigned n;
    pvi_t *pvi;
} rk_t;

/* create 2nd order Runge-Kutta problem */
rk_t *new_rk(pvi_t *pvi, double a, double b, unsigned n);

/* delete 2nd order Runge-Kutta problem */
void free_rk(rk_t *rk);

/* print Runge-Kutta ofof  second order */
void print_rk(rk_t *rk);

/* solve the problem of initial value with Runge-Kutta of some order */
/* return the vector of _y_'s */
vector_t *solve_rk(rk_t *rk, unsigned order);

#endif /* end of include guard: DIFFEQ_H */
