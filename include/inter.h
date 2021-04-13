#ifndef _INTER_H
#define _INTER_H

#include <math.h>

#include "linAlg.h"
#include "matrix.h"
#include "polyn.h"
#include "vector.h"

/* solves minimum square problem with cholesky method */
polyn_t *approx_polyn(vector_t *_x, vector_t *_y, unsigned n, int *stt);

polyn_t *matrix_dif_div_to_polyn(matrix_t *mat);

/* diferencias divididas */
matrix_t *dif_div(matrix_t *pts);

double evaluate_matrix_inter_dif_div(matrix_t *coef, matrix_t *val, double x);

/* cubic splin struct */
typedef struct {
    unsigned n; /* lenght of the splin */
    vector_t *x;
    vector_t *f;
    vector_t *M;
} splin3_t;

/* create splin3_t struct for finding cubic splin to interpolate a function */
/* we asume that the x_values are in the first column, and that are ordered in
 * crecient order, and not containing NaN values */
splin3_t *new_splin3(matrix_t *xf);

/* delete splin3_t struct for finding cubic splin to interpolate a function */
void free_splin3(splin3_t *splin);

/* solve the problem of a cubic splin */
void solve_splin3(splin3_t *splin, double d0, double dn, int *stt);

/* solve the problem of a natural cubic splin */
void solve_nat_splin3(splin3_t *splin, int *stt);

/* evaluate splin in x */
double eval_splin3(splin3_t *splin, double _x);

#endif
