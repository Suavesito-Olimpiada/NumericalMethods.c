#ifndef _VECTOR_H
#define _VECTOR_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vecAux.h"

/* vector type */
typedef struct {
    double *vector;
    unsigned dim;
} vector_t;

/* vector create/free */
vector_t *new_vector(unsigned dim);
void free_vector(vector_t *vec);

/* vector input/output */
vector_t *fread_vector_bin_by_name(char name[]);
vector_t *fread_vector_bin(FILE *pf);
void fwrite_vector_bin_by_name(char name[], vector_t *vec);
void fwrite_vector_bin(FILE *pf, vector_t *vec);

vector_t *fread_vector_text_by_name(char name[]);
vector_t *fread_vector_text(FILE *pf);
void fwrite_vector__text_by_name(char name[], vector_t *vec);
void fwrite_vector_text(FILE *pf, vector_t *vec);

void print_vector(vector_t *vec);
void print_vector_shy(vector_t *vec);

/* vector memory operations */
void vector_copy_i(vector_t *dst, vector_t *src);
vector_t *vector_copy(vector_t *vec);

/* vector fill */
void vector_fill(vector_t *vec, double val);
void vector_fill_zero(vector_t *vec);
void vector_fill_one(vector_t *vec);
void vector_fill_rand(vector_t *vec);

void vector_fill_basis(vector_t *vec, unsigned n);

/* vector operations */
vector_t *vector_transpose(vector_t *vec);

/* vector math operations */
vector_t *vector_sum_vector(vector_t *vec1, vector_t *vec2);
vector_t *vector_sub_vector(vector_t *vec1, vector_t *vec2);
vector_t *vector_mul_vector(vector_t *vec1, vector_t *vec2);
vector_t *vector_div_vector(vector_t *vec1, vector_t *vec2);

vector_t *vector_sum_constant(vector_t *vec, double val);
vector_t *vector_sub_constant(vector_t *vec, double val);
vector_t *vector_scl_constant(vector_t *vec, double val);

void vector_evaluate_function_over_interval(vector_t *vec, double (*f)(double),
                                            double x_0, double x_n);
vector_t *vector_fnc(vector_t *vec, double (*f)(double));

/* in-situ */
void vector_sum_vector_i(vector_t *vec1, vector_t *vec2);
void vector_sub_vector_i(vector_t *vec1, vector_t *vec2);
void vector_mul_vector_i(vector_t *vec1, vector_t *vec2);
void vector_div_vector_i(vector_t *vec1, vector_t *vec2);

void vector_sum_constant_i(vector_t *vec, double val);
void vector_sub_constant_i(vector_t *vec, double val);
void vector_scl_constant_i(vector_t *vec, double val);

void vector_fnc_i(vector_t *vec, double (*f)(double));

/* vector-vector operations */
double vector_dot_product(vector_t *vec1, vector_t *vec2);

/* vector norms */
double vector_l2_norm(vector_t *vec);
double vector_l1_norm(vector_t *vec);
double vector_linf_norm(vector_t *vec);

double vector_lp_norm(vector_t *vec, double p);

#endif
