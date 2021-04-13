#ifndef _LINALG_H
#define _LINALG_H

#include <math.h>
#include <stdbool.h>

#include "matrix.h"
#include "matrix_tridiag.h"
#include "vector.h"

/* linAlg */

extern double tau;

/* in-place copy of the matrix place column */
void vector_col_of_matrix_i(vector_t *vec, matrix_t *mat, unsigned place);
/* return copy of the column number place */
vector_t *vector_col_of_matrix(matrix_t *mat, unsigned place);

/* in-place copy of the matrix diagonal */
void vector_diag_of_matrix_i(vector_t *vec, matrix_t *mat);
/* return copy of the matrix diagonal */
vector_t *vector_diag_of_matrix(matrix_t *mat);

/* matrix times a vector, return a new vector */
vector_t *matrix_mul_vector(matrix_t *mat, vector_t *vec);
vector_t *matrix_tridiag_mul_vector(matrix_tridiag_t *mat, vector_t *vec);

/* matrix multiplication */
matrix_t *matrix_mul_matrix(matrix_t *mat1, matrix_t *mat2);
/* matrix_tridiag_t *matrix_tridiag_mul_matrix_tridiag(matrix_tridiag_t *mat1,
 */
/*                                                     matrix_tridiag_t *mat2);
 */
/* matrix_tridiag_t *matrix_tridiag_mul_matrix(matrix_tridiag_t *mat1, */
/*                                             matrix_t *mat2); */
/* matrix times matrix transpose */
matrix_t *matrix_mul_matrix_t(matrix_t *mat);
/* matrix transpose times matrix */
matrix_t *matrix_t_mul_matrix(matrix_t *mat);
/* matrix times diiaigonal matrix */
bool matrix_mul_matrix_diag_i(matrix_t *mat, vector_t *diag);
matrix_t *matrix_mul_matrix_diag(matrix_t *mat, vector_t *diag);
/* matrix times Givens matrix */
bool matrix_mul_matrix_givens_i(matrix_t *mat, unsigned i, unsigned j,
                           double s, double c);
matrix_t *matrix_mul_matrix_givens(matrix_t *mat, unsigned _i, unsigned _j,
                              double s, double c);

/* error of a linear equation system */
double lin_system_error(matrix_t *A_m, vector_t *b_v, vector_t *x_v);
double lin_system_error_matrix_tridiag(matrix_tridiag_t *mat, vector_t *_b,
                                       vector_t *_x);

/* print error messages corresponding to the error codes */
void lin_print_errors(int stt);

/* solve upper matrix with forward sustitution */
vector_t *lin_up_matrix_solve(matrix_t *mat, vector_t *vec, int *stt);

/* solve upper matrix with backward sustitution */
vector_t *lin_low_matrix_solve(matrix_t *mat, vector_t *vec, int *stt);

/* solve upper matrix with forward sustitution with 1's in the diagonal */
vector_t *lin_up_unit_matrix_solve(matrix_t *mat, vector_t *vec, int *stt);

/* solve lower matrix with backward sustitution with 1's in the diagonal */
vector_t *lin_low_unit_matrix_solve(matrix_t *mat, vector_t *vec, int *stt);

/* using Doolitle finds the LU factorization of mat */
matrix_t *lin_matrix_lu_doolitle(matrix_t *mat, int *stt);

/* solve the system using L*U*x = y */
vector_t *solve_lu_doolitle(matrix_t *mat, vector_t *vec, int *stt);

/* solution struct for LU with pivoting */
typedef struct {
    matrix_t *L;
    matrix_t *U;
    unsigned *p;
} lup_t;

/* create lup_t struct for solving with LU with pivoting */
lup_t *new_lup(matrix_t *mat);

/* delete lup_t struct for solving with LU with pivoting */
void free_lup(lup_t *lup);

/* find the LU factorization with partial pivoting of mat */
lup_t *lin_matrix_lu_pivote(matrix_t *mat, int *stt);

/* generate the solution of a system using L*U*x = P*y */
vector_t *gen_solve_lup(lup_t *lup, vector_t *vec, int *stt);

/* solve the system using L*U*x = P*y */
vector_t *solve_lup(matrix_t *mat, vector_t *vec, int *stt);

/* find the Cholesky factorization of mat */
matrix_t *lin_matrix_cholesky(matrix_t *mat, int *stt);

/* solve the system using L*L^t*x = y */
vector_t *solve_cholesky(matrix_t *mat, vector_t *vec, int *stt);

/* tridiag matrix solving method */
vector_t *solve_matrix_tridiag(matrix_tridiag_t *mat, vector_t *vec, int *stt);

/* tridiag matrix solving with the Jacobi method */
vector_t *solve_jacobi_matrix_tridiag(matrix_tridiag_t *mat, vector_t *vec,
                                      unsigned max_iter, int *stt);

/* solution struct for Gauss-Seidel */
typedef struct {
    vector_t *x;
    double error;
    unsigned iter;
} gs_t;

/* create gs_t struct for solving with Gauss-Seidel */
gs_t *new_gs(unsigned dim);

/* delete gs_t struct for solving with Gauss-Seidel */
void free_gs(gs_t *gs);

/* print values of the Gauss-Seidel solution */
void print_gs(gs_t *gs);

/* Single in-place step in Gauss Seidel method for tridiagonal matrix */
void lin_gauss_seidel_matrix_tridiag_i(matrix_tridiag_t *mat, vector_t *_b,
                                       vector_t *_x, int *stt);

/* tridiag matrix solving with the Jacobi method */
gs_t *solve_gauss_seidel_matrix_tridiag(matrix_tridiag_t *mat, vector_t *vec,
                                        unsigned max_iter, int *stt);

/* solves minimum square error problem with cholesky method */
vector_t *solve_mse_cholesky(matrix_t *mat, vector_t *vec, int *stt);

/* solution struct for Jacobi iterative */
typedef struct {
    matrix_t *V;
    vector_t *d;
    double bmax;
    unsigned iter;
} ji_t;

/* create ji_t struct for Jacobi matrix eigen-value decomposition */
ji_t *new_ji(unsigned dim);

/* delete ji_t struct for Jacobi matrix eigen-value decomposition */
void free_ji(ji_t *ji);

/* print values of the Jacobi matrix eigen-value decomposition */
void print_ji(ji_t *ji);

/* matrix iterative Jacobi method for eigen-value decomposition */
ji_t *iterative_jacobi_matrix_eigen(matrix_tridiag_t *mat, unsigned max_iter,
                                    int *stt);

#endif
