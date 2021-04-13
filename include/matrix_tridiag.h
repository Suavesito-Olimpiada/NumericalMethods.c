#ifndef _MATRIX_TRIDIAG_H
#define _MATRIX_TRIDIAG_H 1

#include "matrix.h"

/* matrix-tridiag type */
typedef matrix_t matrix_tridiag_t;

/* matrix-tridiag create/free */
matrix_tridiag_t *new_matrix_tridiag(unsigned dim);
void free_matrix_tridiag(matrix_tridiag_t *mat);

/* matrix-tridiag input/output */
matrix_tridiag_t *fread_matrix_tridiag_bin_by_name(char *name);
matrix_tridiag_t *fread_matrix_tridiag_bin(FILE *pf);
void fwrite_matrix_tridiag_bin_by_name(char *name, matrix_tridiag_t *mat);
void fwrite_matrix_tridiag_bin(FILE *pf, matrix_tridiag_t *mat);

matrix_tridiag_t *fread_matrix_tridiag_text_by_name(char *name);
matrix_tridiag_t *fread_matrix_tridiag_text(FILE *pf);
void fwrite_matrix_tridiag_text_by_name(char *name, matrix_tridiag_t *mat);
void fwrite_matrix_tridiag_text(FILE *pf, matrix_tridiag_t *mat);

void print_matrix_tridiag(matrix_tridiag_t *mat);
void print_matrix_tridiag_shy(matrix_tridiag_t *mat);

/* matrix-tridiag operations */
matrix_tridiag_t *matrix_tridiag_copy(matrix_tridiag_t *mat);
matrix_t *matrix_tridiag_to_matrix(matrix_tridiag_t *mat);
matrix_tridiag_t *matrix_to_matrix_tridiag(matrix_t *mat);

/* matrix-tridiag fill */
void matrix_tridiag_fill(matrix_tridiag_t *mat, double val);
void matrix_tridiag_fill_zero(matrix_tridiag_t *mat);
void matrix_tridiag_fill_rand(matrix_tridiag_t *mat);

void matrix_tridiag_fill_rand_up(matrix_tridiag_t *mat);
void matrix_tridiag_fill_rand_low(matrix_tridiag_t *mat);

/* matrix-tridiag operations */
void matrix_tridiag_transpose_i(matrix_tridiag_t *mat);
matrix_tridiag_t *matrix_tridiag_transpose(matrix_tridiag_t *mat);

/* check if the matrix mat is dominant diagonal */
bool matrix_tridiag_is_diagonal_dominant(matrix_tridiag_t *mat);

/* matrix-tridiag math operations */
matrix_tridiag_t *matrix_tridiag_sum_matrix_tridiag(matrix_tridiag_t *mat1,
                                                    matrix_tridiag_t *mat2);
matrix_tridiag_t *matrix_tridiag_sub_matrix_tridiag(matrix_tridiag_t *mat1,
                                                    matrix_tridiag_t *mat2);
matrix_tridiag_t *matrix_tridiag_div_matrix_tridiag(matrix_tridiag_t *mat1,
                                                    matrix_tridiag_t *mat2);

matrix_t *matrix_tridiag_sum_matrix(matrix_tridiag_t *mat1, matrix_t *mat2);
matrix_t *matrix_tridiag_sub_matrix(matrix_tridiag_t *mat1, matrix_t *mat2);
matrix_tridiag_t *matrix_tridiag_div_matrix(matrix_tridiag_t *mat1,
                                            matrix_t *mat2);

matrix_t *matrix_sub_matrix_tridiag(matrix_t *mat1, matrix_tridiag_t *mat2);

matrix_tridiag_t *matrix_tridiag_sum_constant(matrix_tridiag_t *mat,
                                              double val);
matrix_tridiag_t *matrix_tridiag_scl_constant(matrix_tridiag_t *mat,
                                              double val);

/* matrix_tridiag norms */
double matrix_tridiag_nat1_norm(matrix_tridiag_t *mat);
double matrix_tridiag_natinf_norm(matrix_tridiag_t *mat);

double matrix_tridiag_frobenius_norm(matrix_tridiag_t *mat);
double matrix_tridiag_max_norm(matrix_tridiag_t *mat);

#endif
