#ifndef _MATRIX_H
#define _MATRIX_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vecAux.h"

/* matrix type */
typedef struct {
    double *data;
    double **matrix;
    unsigned dim1, dim2;
} matrix_t;

/* matrix create/free */
matrix_t *new_matrix(unsigned dim1, unsigned dim2);
void free_matrix(matrix_t *mat);

/* matrix input/output */
matrix_t *fread_matrix_bin_by_name(char *name);
matrix_t *fread_matrix_bin(FILE *pf);
void fwrite_matrix_bin_by_name(char *name, matrix_t *mat);
void fwrite_matrix_bin(FILE *pf, matrix_t *mat);

matrix_t *fread_matrix_text_by_name(char *name);
matrix_t *fread_matrix_text(FILE *pf);
void fwrite_matrix_text_by_name(char *name, matrix_t *mat);
void fwrite_matrix_text(FILE *pf, matrix_t *mat);

void print_matrix(matrix_t *mat);
void print_matrix_shy(matrix_t *mat);

/* matrix operations */
void matrix_copy_i(matrix_t *dst, matrix_t *src);
matrix_t *matrix_copy(matrix_t *mat);

/* matrix fill */
void matrix_fill(matrix_t *mat, double val);
void matrix_fill_zero(matrix_t *mat);
void matrix_fill_rand(matrix_t *mat);

void matrix_fill_rand_up(matrix_t *mat);
void matrix_fill_rand_low(matrix_t *mat);

void matrix_fill_identity(matrix_t *mat);
void matrix_fill_rand_diag(matrix_t *mat);

/* matrix operations */
void matrix_transpose_i(matrix_t *mat);
matrix_t *matrix_transpose(matrix_t *mat);

/* check iff the matrix mat is symmetric with tolerance _tau */
bool matrix_is_sym(matrix_t *mat, double _tau);

/* matrix math operations */
matrix_t *matrix_sum_matrix(matrix_t *mat1, matrix_t *mat2);
matrix_t *matrix_sub_matrix(matrix_t *mat1, matrix_t *mat2);
matrix_t *matrix_div_matrix(matrix_t *mat1, matrix_t *mat2);

matrix_t *matrix_sum_constant(matrix_t *mat, double val);
matrix_t *matrix_sub_constant(matrix_t *mat, double val);
matrix_t *matrix_scl_constant(matrix_t *mat, double val);

matrix_t *matrix_fnc(matrix_t *mat, double (*f)(double));

/* in-situ */
void matrix_sum_matrix_i(matrix_t *mat1, matrix_t *mat2);
void matrix_sub_matrix_i(matrix_t *mat1, matrix_t *mat2);
void matrix_div_matrix_i(matrix_t *mat1, matrix_t *mat2);

void matrix_sum_constant_i(matrix_t *mat, double val);
void matrix_sub_constant_i(matrix_t *mat, double val);
void matrix_scl_constant_i(matrix_t *mat, double val);

void matrix_fnc_i(matrix_t *mat, double (*f)(double));

/* matrix norms */
double matrix_nat1_norm(matrix_t *mat);
double matrix_natinf_norm(matrix_t *mat);

double matrix_frobenius_norm(matrix_t *mat);
double matrix_max_norm(matrix_t *mat);

/* matrix utils */
/* find the biggest element above the diagonal */
double matrix_max_diag(matrix_t *mat, unsigned *i, unsigned *j);

#endif
