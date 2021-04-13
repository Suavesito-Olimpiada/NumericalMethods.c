#include "linAlg.h"

void vector_col_of_matrix_i(vector_t *vec, matrix_t *mat, unsigned place)
{
    if (!mat || !vec || vec->dim != mat->dim2 || place > mat->dim2)
        return;

    unsigned i = 0;
    for (i = 0; i < vec->dim; ++i) vec->vector[i] = mat->matrix[i][place];
}

vector_t *vector_col_of_matrix(matrix_t *mat, unsigned place)
{
    if (!mat || place > mat->dim2)
        return NULL;

    vector_t *new = new_vector(mat->dim1);
    if (!new)
        return NULL;

    vector_col_of_matrix_i(new, mat, place);

    return new;
}

void vector_diag_of_matrix_i(vector_t *vec, matrix_t *mat)
{
    if (!mat || !vec)
        return;

    if (vec->dim != (mat->dim1 < mat->dim2 ? mat->dim1 : mat->dim2))
        return;

    unsigned i = 0;
    for (i = 0; i < vec->dim; ++i) {
        vec->vector[i] = mat->matrix[i][i];
    }
}

vector_t *vector_diag_of_matrix(matrix_t *mat)
{
    if (!mat)
        return NULL;

    vector_t *new = new_vector((mat->dim1 < mat->dim2) ? mat->dim1 : mat->dim2);
    if (!new)
        return NULL;

    vector_diag_of_matrix_i(new, mat);

    return new;
}

vector_t *matrix_mul_vector(matrix_t *mat, vector_t *vec)
{
    if (!mat || !vec)
        return NULL;

    if (mat->dim2 != vec->dim)
        return NULL;

    vector_t *new = new_vector(mat->dim1);
    if (!new)
        return NULL;

    unsigned i = 0;
    double *mat_data = mat->data, *vec_data = vec->vector;
    for (i = 0; i < mat->dim1; ++i, mat_data += mat->dim2)
        new->vector[i] = __dot_product(mat_data, vec_data, vec->dim);

    return new;
}

vector_t *matrix_tridiag_mul_vector(matrix_tridiag_t *mat, vector_t *vec)
{
    if (!mat || !vec)
        return NULL;

    if (mat->dim1 != vec->dim)
        return NULL;

    vector_t *new = new_vector(mat->dim1);
    if (!new)
        return NULL;

    unsigned n = mat->dim1;
    unsigned i = 0;
    double **A = mat->matrix;
    double *x = vec->vector;
    double *b = new->vector;
    b[0] = (A[0][1] * x[0]) + (A[0][2] * x[1]);
    for (i = 1; i < n - 1; ++i) {
        b[i] = (A[i][0] * x[i - 1]) + (A[i][1] * x[i]) + (A[i][2] * x[i + 1]);
    }
    b[n - 1] = (A[n - 1][0] * x[n - 2]) + (A[n - 1][1] * x[n - 1]);

    return new;
}

matrix_t *matrix_mul_matrix(matrix_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if ((mat1->dim1 != mat2->dim1) || (mat1->dim2 != mat2->dim2))
        return NULL;

    matrix_t *new = new_matrix(mat1->dim1, mat2->dim2);
    if (!new)
        return NULL;

    matrix_t *trans = matrix_transpose(mat2);

    double **A = new->matrix;
    unsigned i = 0, j = 0;
    for (i = 0; i < mat1->dim1; ++i) {
        for (j = 0; j < mat2->dim2; ++j) {
            A[i][j] =
                __dot_product(mat1->data + (i * mat1->dim2),
                              trans->data + (j * trans->dim2), mat1->dim2);
        }
    }

    free_matrix(trans);

    return new;
}

/* matrix_tridiag_t *matrix_tridiag_mul_matrix_tridiag(matrix_tridiag_t *mat1,
                                                    matrix_tridiag_t *mat2)
{
    if(!mat1 || !mat2)
        return NULL;

    if(mat->dim1 != mat2->dim1)
        return NULL;



    return matrix_mul_matrix((matrix_t *)mat1, (matrix_t *)mat2);
}

matrix_tridiag_t *matrix_tridiag_mul_matrix(matrix_tridiag_t *mat1,
                                            matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if (mat1->dim1 != mat2->dim1 || mat1->dim1 != mat2->dim2)
        return NULL;

    matrix_t *new = matrix_to_matrix_tridiag(mat2);
    if (!new)
        return NULL;

    unsigned i = 0;
    double *it = mat2->data, *it2 = new->data;
    for (i = 1; i < mat2->dim1 * 3 - 1; ++i) *it2 *= *it;

    return new;
} */

matrix_t *matrix_mul_matrix_t(matrix_t *mat)
{
    if (!mat)
        return NULL;

    matrix_t *new = new_matrix(mat->dim1, mat->dim1);
    if (!new)
        return NULL;

    double **M = new->matrix;
    unsigned i = 0, j = 0;
    for (i = 0; i < mat->dim1; ++i) {
        for (j = 0; j < mat->dim1; ++j) {
            M[i][j] = __dot_product(mat->data + (i * mat->dim2),
                                    mat->data + (j * mat->dim2), mat->dim2);
        }
    }

    return new;
}

matrix_t *matrix_t_mul_matrix(matrix_t *mat)
{
    if (!mat)
        return NULL;

    /* first we transpose the matrix, so the calculation can be made row-major
     * and take advantange of cpu cache. It can be faster for large matrix to
     * take this advantange over not having to copy the matrix */
    matrix_t *trans = matrix_transpose(mat);
    if (!trans)
        return NULL;

    matrix_t *new = matrix_mul_matrix_t(trans);

    /* free transposed matrix */
    free_matrix(trans);

    return new;
}

bool matrix_mul_matrix_diag_i(matrix_t *mat, vector_t *diag)
{
    if (!mat)
        return false;

    if (!diag)
        return false;

    if (mat->dim1 != mat->dim2 || mat->dim1 != diag->dim)
        return false;

    unsigned n = mat->dim1;
    unsigned i = 0, j = 0;
    double *d = diag->vector;
    double **A = mat->matrix;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            A[i][j] *= d[j];
        }
    }

    return true;
}

matrix_t *matrix_mul_matrix_diag(matrix_t *mat, vector_t *diag)
{
    if (!mat)
        return NULL;

    if (!diag)
        return NULL;

    if (mat->dim1 != mat->dim2 || mat->dim1 != diag->dim)
        return NULL;

    matrix_t *new = matrix_copy(mat);
    if (!new)
        return NULL;

    matrix_mul_matrix_diag_i(new, diag);

    return new;
}

bool matrix_mul_matrix_givens_i(matrix_t *mat, unsigned i, unsigned j,
                                double s, double c)
{
    if (!mat)
        return false;

    if (mat->dim2 <= i || mat->dim2 <= j || i == j)
        return false;

    double a_i = 0.0;
    double a_j = 0.0;
    unsigned k = 0;
    double **A = mat->matrix;
    for (k = 0; k < mat->dim1; ++k) {
        a_i = A[k][i];
        a_j = A[k][j];
        A[k][i] = c * a_i - s * a_j;
        A[k][j] = s * a_i + c * a_j;
    }

    return true;
}

matrix_t *matrix_mul_matrix_givens(matrix_t *mat, unsigned _i, unsigned _j,
                                   double s, double c)
{
    if (!mat)
        return NULL;

    if (mat->dim2 <= _i || mat->dim2 <= _j || _i == _j)
        return NULL;

    if (_i > _j) {
        unsigned temp = _i;
        _i = _j;
        _j = temp;
    }

    matrix_t *new = matrix_copy(mat);
    if (!new)
        return NULL;

    matrix_mul_matrix_givens_i(new, _i, _j, s, c);

    return new;
}

double lin_system_error(matrix_t *A_m, vector_t *b_v, vector_t *x_v)
{
    if (!A_m || !b_v || !x_v)
        return NAN;

    if ((A_m->dim1 != b_v->dim) || (A_m->dim2 != x_v->dim))
        return NAN;

    double er = 0.0;

    unsigned i = 0;
    double **it = A_m->matrix, *it2 = b_v->vector;
    for (i = 0; i < A_m->dim1; ++i, ++it, ++it2)
        er += pow(__dot_product(*it, x_v->vector, A_m->dim2) - *it2, 2.0);

    er = sqrt(er);

    return er;
}

double lin_system_error_matrix_tridiag(matrix_tridiag_t *mat, vector_t *_b,
                                       vector_t *_x)
{
    if (!mat || !_b || !_x)
        return NAN;

    if ((mat->dim1 != _b->dim) || (mat->dim1 != _x->dim))
        return NAN;

    double error = 0.;
    unsigned n = mat->dim1;
    unsigned i = 0;
    double **A = mat->matrix;
    double *b = _b->vector;
    double *x = _x->vector;
    error += pow(b[0] - ((A[0][1] * x[0]) + (A[0][2] * x[1])), 2.);
    for (i = 1; i < n - 1; ++i) {
        error += pow(b[i] - ((A[i][0] * x[i - 1]) + (A[i][1] * x[i]) +
                             (A[i][2] * x[i + 1])),
                     2.);
    }
    error += pow(
        b[n - 1] - ((A[n - 1][0] * x[n - 2]) + (A[n - 1][1] * x[n - 1])), 2.);

    error = sqrt(error);

    return error;
}

void lin_print_errors(int stt)
{
    switch (stt) {
        case -10:
            perror("Memory error: NULL argument.\n");
            break;
        case -11:
            perror("Memory error: returning vector.\n");
            break;
        case -12:
            perror("Memory error: returning matrix.\n");
            break;
        case -13:
            perror("Memory error: internal variables.\n");
            break;
        case 11:
            perror("Numerical error: incompatible dimentions.\n");
            break;
        case 12:
            perror("Numerical error: incompatible matrix.\n");
            break;
        case 13:
            perror("Numerical error: dividing by zero (matrix value).\n");
            break;
        case 14:
            perror("Numerical error: NaN or Inf values.\n");
            break;
        case 15:
            perror("Numerical error: Maximum number of iterations.\n");
            break;
    }
}

vector_t *lin_up_matrix_solve(matrix_t *mat, vector_t *vec, int *stt)
{
    *stt = -10; /* memory failure (-1) */
                /* no memory for original object (0) */
    if (!mat || !vec || !stt)
        return NULL;

    *stt = 11; /* no-solution (1), incopatible dimentions (1) */
    if (mat->dim1 != vec->dim)
        return NULL;

    *stt = -11; /* memory failure (-1) */
                /* no memory for solution vector (1) */
    vector_t *new = new_vector(mat->dim2);
    if (!new)
        return NULL;

    *stt = 0; /* everything ok */

    if (mat->dim2 > mat->dim1)
        *stt = 2; /* multiple solutions (2) */

    if (mat->dim1 > mat->dim2) {
        *stt = 12; /* no-solution (1), incopatible numbers (2) */
        unsigned i = 0;
        for (i = vec->dim - 1; i > mat->dim1 - mat->dim2; --i) {
            if (vec->vector[i] != 0) {
                free_vector(new);
                return NULL;
            }
        }
    }

    vector_fill_zero(new);

    double **A = mat->matrix, *b = vec->vector, *x = new->vector;
    double v;
    unsigned i = 0;
    for (i = 0; i < vec->dim; ++i) {
        if (fabs(A[vec->dim - i - 1][vec->dim - i - 1]) < tau) {
            *stt = 13; /* no-solution (1), zero in diagonal (3) */
            free_vector(new);
            return NULL;
        }
        v = __dot_product(A[vec->dim - i - 1], x, vec->dim);
        x[vec->dim - i - 1] =
            (b[vec->dim - i - 1] - v) / A[vec->dim - i - 1][vec->dim - i - 1];
        if (isinf(v) || isnan(v) || isinf(x[vec->dim - i - 1]) ||
            isnan(x[vec->dim - i - 1])) {
            *stt = 14; /* no-solution (1), numerical error (4) */
            free_vector(new);
            return NULL;
        }
    }

    return new;
}

vector_t *lin_low_matrix_solve(matrix_t *mat, vector_t *vec, int *stt)
{
    *stt = -10; /* memory failure (-1) */
                /* no memory for original object (0) */
    if (!mat || !vec || !stt)
        return NULL;

    *stt = 11; /* no-solution (1), incopatible dimentions (1) */
    if (mat->dim1 != vec->dim)
        return NULL;

    *stt = -11; /* memory failure (-1) */
                /* no memory for solution vector (1) */
    vector_t *new = new_vector(mat->dim2);
    if (!new)
        return NULL;

    *stt = 0; /* everything ok */

    if (mat->dim2 > mat->dim1)
        *stt = 2; /* multiple solutions (2) */

    vector_fill_zero(new);

    double **A = mat->matrix, *b = vec->vector, *x = new->vector;
    double v;
    unsigned i = 0;
    for (i = 0; i < vec->dim; ++i) {
        if (fabs(A[i][i]) < tau) {
            *stt = 13; /* no-solution (1), zero in diagonal (3) */
            free_vector(new);
            return NULL;
        }
        v = __dot_product(A[i], x, vec->dim);
        x[i] = (b[i] - v) / A[i][i];
        if (isinf(v) || isnan(v) || isinf(x[i]) || isnan(x[i])) {
            *stt = 14; /* no-solution (1), numerical error (4) */
            free_vector(new);
            return NULL;
        }
    }

    return new;
}

vector_t *lin_up_unit_matrix_solve(matrix_t *mat, vector_t *vec, int *stt)
{
    *stt = -10; /* memory failure (-1) */
                /* no memory for original object (0) */
    if (!mat || !vec || !stt)
        return NULL;

    *stt = 11; /* no-solution (1), incopatible dimentions (1) */
    if (mat->dim1 != vec->dim)
        return NULL;

    *stt = -11; /* memory failure (-1) */
                /* no memory for solution vector (1) */
    vector_t *new = new_vector(mat->dim2);
    if (!new)
        return NULL;

    *stt = 0; /* everything ok */

    if (mat->dim2 > mat->dim1)
        *stt = 2; /* multiple solutions (2) */

    if (mat->dim1 > mat->dim2) {
        *stt = 12; /* no-solution (1), incopatible numbers (2) */
        unsigned i = 0;
        for (i = vec->dim - 1; i > mat->dim1 - mat->dim2; --i) {
            if (vec->vector[i] != 0) {
                free_vector(new);
                return NULL;
            }
        }
    }

    vector_fill_zero(new);

    double **A = mat->matrix, *b = vec->vector, *x = new->vector;
    double v;
    unsigned i = 0;
    for (i = 0; i < vec->dim; ++i) {
        v = __dot_product(A[vec->dim - i - 1] + i, x + i, vec->dim - i);
        x[vec->dim - i - 1] = b[vec->dim - i - 1] - v;
        if (isinf(v) || isnan(v) || isinf(x[vec->dim - i - 1]) ||
            isnan(x[vec->dim - i - 1])) {
            *stt = 14; /* no-solution (1), numerical error (4) */
            free_vector(new);
            return NULL;
        }
    }

    return new;
}

vector_t *lin_low_unit_matrix_solve(matrix_t *mat, vector_t *vec, int *stt)
{
    *stt = -10; /* memory failure (-1) */
                /* no memory for original object (0) */
    if (!mat || !vec || !stt)
        return NULL;

    *stt = 11; /* no-solution (1), incopatible dimentions (1) */
    if (mat->dim1 != vec->dim)
        return NULL;

    *stt = -11; /* memory failure (-1) */
                /* no memory for solution vector (1) */
    vector_t *new = new_vector(mat->dim2);
    if (!new)
        return NULL;

    *stt = 0; /* everything ok */

    if (mat->dim2 > mat->dim1)
        *stt = 2; /* multiple solutions (2) */

    vector_fill_zero(new);

    double **A = mat->matrix, *b = vec->vector, *x = new->vector;
    double v;
    unsigned i = 0;
    for (i = 0; i < vec->dim; ++i) {
        if (fabs(A[i][i]) < tau) {
            *stt = 13; /* no-solution (1), zero in diagonal (3) */
            free_vector(new);
            return NULL;
        }
        v = __dot_product(A[i], x, i) + x[i];
        x[i] = b[i] - v;
        if (isinf(v) || isnan(v) || isinf(x[i]) || isnan(x[i])) {
            *stt = 14; /* no-solution (1), numerical error (4) */
            free_vector(new);
            return NULL;
        }
    }

    return new;
}

matrix_t *lin_matrix_lu_doolitle(matrix_t *mat, int *stt)
{
    *stt = -10;
    if (!mat)
        return NULL;

    *stt = 11;
    if (mat->dim1 != mat->dim2)
        return NULL;

    *stt = -12;
    matrix_t *new = matrix_copy(mat);
    if (!new)
        return NULL;

    unsigned i = 0, j = 0, k = 0;
    double sum = 0.0;
    double **ALU = new->matrix;
    for (j = 0; j < mat->dim2; ++j) {
        for (i = 0; i <= j; ++i) {
            sum = 0.0;
            for (k = 0; k < i; ++k) sum += ALU[i][k] * ALU[k][j];
            ALU[i][j] = ALU[i][j] - sum;
            if (isinf(ALU[i][j]) || isnan(ALU[i][j])) {
                *stt = 14;
                free_matrix(new);
                return NULL;
            }
        }
        if (fabs(ALU[j][j]) < tau) {
            *stt = 13;
            free_matrix(new);
            return NULL;
        }
        for (i = j + 1; i < mat->dim1; ++i) {
            sum = 0.0;
            for (k = 0; k < j; ++k) sum += ALU[i][k] * ALU[k][j];
            ALU[i][j] = (1.0 / ALU[j][j]) * (ALU[i][j] - sum);
        }
    }
    *stt = 0;

    return new;
}

vector_t *solve_lu_doolitle(matrix_t *mat, vector_t *vec, int *stt)
{
    matrix_t *ALU = lin_matrix_lu_doolitle(mat, stt);
    if (!ALU) {
        return NULL;
    }

    /* solve the system L * y = vec */
    vector_t *y = lin_low_unit_matrix_solve(ALU, vec, stt);
    if (!y) {
        free_matrix(ALU);
        return NULL;
    }

    /* solve the system U * x = y */
    vector_t *x = lin_up_matrix_solve(ALU, y, stt);
    if (!x) {
        free_vector(y);
        free_matrix(ALU);
        return NULL;
    }

    free_vector(y);
    free_matrix(ALU);

    return x;
}

lup_t *new_lup(matrix_t *mat)
{
    if (!mat)
        return NULL;

    if (mat->dim1 != mat->dim2)
        return NULL;

    lup_t *new = (lup_t *)malloc(sizeof(lup_t));
    if (!new)
        goto f3;

    new->L = new_matrix(mat->dim1, mat->dim2);
    if (!new->L)
        goto f2;

    new->U = new_matrix(mat->dim1, mat->dim2);
    if (!new->U)
        goto f1;

    new->p = (unsigned *)malloc(mat->dim1 * sizeof(unsigned));
    if (!new->p)
        goto f;

    matrix_copy_i(new->U, mat);
    matrix_fill_identity(new->L);

    unsigned i = 0;
    for (i = 0; i < mat->dim1; ++i) new->p[i] = i;

    return new;

f:
    free_matrix(new->U);
f1:
    free_matrix(new->L);
f2:
    free_lup(new);
f3:
    return new;
}

void free_lup(lup_t *lup)
{
    if (!lup)
        return;
    if (lup->L)
        free_matrix(lup->L);
    if (lup->U)
        free_matrix(lup->U);
    if (lup->p)
        free(lup->p);
    free(lup);
}

lup_t *lin_matrix_lu_pivote(matrix_t *mat, int *stt)
{
    *stt = -10;
    if (!mat)
        return NULL;

    *stt = 11;
    if (mat->dim1 != mat->dim2)
        return NULL;

    *stt = -12;
    lup_t *new = new_lup(mat);
    if (!new)
        return NULL;

    unsigned i = 0, j = 0, k = 0;
    unsigned n = mat->dim1;
    unsigned r = 0;
    unsigned *p = new->p;
    double **L = new->L->matrix;
    double **U = new->U->matrix;
    double faux;
    unsigned aux;
    for (k = 0; k < n; ++k) {
        r = k;
        faux = fabs(U[r][k]);
        for (i = 0; i < n; ++i) {
            r = fabs(U[r][k]) > faux ? i : r;
        }

        if (fabs(U[r][k]) < tau) {
            *stt = 13;
            free_lup(new);
            return NULL;
        }

        if (r != k) {
            for (i = 0; i < n; ++i) {
                faux = U[k][i];
                U[k][i] = U[r][i];
                U[r][i] = faux;
            }

            aux = p[k];
            p[k] = p[r];
            p[r] = aux;

            if (k > 0) {
                for (i = 0; i < k; ++i) {
                    faux = L[k][i];
                    L[k][i] = L[r][i];
                    L[r][i] = faux;
                }
            }
        }

        for (i = k + 1; i < n; ++i) {
            L[i][k] = U[i][k] / U[k][k];
            for (j = k; j < n; ++j) {
                U[i][j] = U[i][j] - L[i][k] * U[k][j];
            }
        }
    }

    *stt = 0;

    return new;
}

vector_t *gen_solve_lup(lup_t *lup, vector_t *vec, int *stt)
{
    *stt = -11;
    vector_t *_vec = vector_copy(vec);
    if (!_vec)
        return NULL;

    unsigned i = 0;
    for (i = 0; i < _vec->dim; ++i) {
        _vec->vector[i] = vec->vector[lup->p[i]];
    }

    /* solve the system L * y = vec */
    vector_t *y = lin_low_unit_matrix_solve(lup->L, _vec, stt);

    if (!y) {
        free_vector(_vec);
        return NULL;
    }

    /* solve the system U * x = y */
    vector_t *x = lin_up_matrix_solve(lup->U, y, stt);

    free_vector(_vec);
    free_vector(y);

    return x;
}

vector_t *solve_lup(matrix_t *mat, vector_t *vec, int *stt)
{
    lup_t *lup = lin_matrix_lu_pivote(mat, stt);
    if (!lup) {
        return NULL;
    }

    vector_t *x = gen_solve_lup(lup, vec, stt);

    free_lup(lup);

    return x;
}

matrix_t *lin_matrix_cholesky(matrix_t *mat, int *stt)
{
    *stt = -10;
    if (!mat)
        return NULL;

    *stt = 11;
    if (mat->dim1 != mat->dim2)
        return NULL;

    *stt = -12;
    matrix_t *new = new_matrix(mat->dim1, mat->dim2);
    if (!new)
        return NULL;
    matrix_fill_zero(new);

    unsigned i = 0, j = 0, k = 0;
    unsigned n = mat->dim1;
    double **A = mat->matrix;
    double **L = new->matrix;
    double faux;
    for (j = 0; j < n; ++j) {
        faux = 0.0;
        for (k = 0; k < j; ++k) {
            faux += L[j][k] * L[j][k];
        }

        faux = A[j][j] - faux;
        if (faux < 0.0) {
            *stt = 14;
            free_matrix(new);
            return NULL;
        }

        L[j][j] = sqrt(faux);
        if (fabs(L[j][j]) < tau) {
            *stt = 13;
            free_matrix(new);
            return NULL;
        }

        for (i = j + 1; i < n; ++i) {
            faux = 0.0;
            for (k = 0; k < j; ++k) {
                faux += L[i][k] * L[j][k];
            }
            L[i][j] = (1.0 / L[j][j]) * (A[i][j] - faux);
        }
    }

    *stt = 0;

    return new;
}

vector_t *solve_cholesky(matrix_t *mat, vector_t *vec, int *stt)
{
    if (!mat || !vec)
        return NULL;

    *stt = 12;
    if (!matrix_is_sym(mat, tau))
        return NULL;

    matrix_t *chol = lin_matrix_cholesky(mat, stt);
    if (!chol) {
        return NULL;
    }

    /* solve the system L * y = vec */
    vector_t *y = lin_low_matrix_solve(chol, vec, stt);

    if (!y) {
        free_matrix(chol);
        return NULL;
    }

    matrix_transpose_i(chol);

    /* solve the system U * x = y */
    vector_t *x = lin_up_matrix_solve(chol, y, stt);

    free_matrix(chol);
    free_vector(y);

    return x;
}

vector_t *solve_matrix_tridiag(matrix_tridiag_t *mat, vector_t *vec, int *stt)
{
    *stt = -10;
    if (!mat || !vec)
        return NULL;

    *stt = 11;
    if (mat->dim1 != vec->dim)
        return NULL;

    *stt = 11;
    if (mat->dim2 != 3)
        return NULL;

    *stt = -11;
    vector_t *new = new_vector(vec->dim);
    if (!new)
        return NULL;

    *stt = -13;

    double *b = (double *)malloc(vec->dim * sizeof(double));
    if (b == NULL) {
        free_vector(new);
        return NULL;
    }

    double *c = (double *)malloc(vec->dim * sizeof(double));
    if (c == NULL) {
        free(b);
        free_vector(new);
        return NULL;
    }

    double *d = (double *)malloc(vec->dim * sizeof(double));
    if (d == NULL) {
        free(c);
        free(b);
        free_vector(new);
        return NULL;
    }

    *stt = 0;

    double **A = mat->matrix, *v = vec->vector, *x = new->vector;
    memset(b, 0, vec->dim * sizeof(double));
    memset(c, 0, vec->dim * sizeof(double));
    memset(d, 0, vec->dim * sizeof(double));
    b[0] = A[0][1], c[0] = A[0][2], d[0] = v[0];
    unsigned i = 0;
    for (i = 1; i < vec->dim; ++i) {
        b[i] = b[i - 1] * A[i][1] - A[i][0] * c[i - 1];
        c[i] = b[i - 1] * A[i][2];
        d[i] = b[i - 1] * v[i] - A[i][0] * d[i - 1];
    }
    if (fabs(b[vec->dim - 1]) < tau || isnan(b[vec->dim - 1])) {
        *stt = isnan(b[vec->dim - 1]) ? 14 : 13;

        free_vector(new);
        free(b), free(c), free(d);

        return NULL;
    }
    x[vec->dim - 1] = d[vec->dim - 1] / b[vec->dim - 1];
    for (i = 1; i < vec->dim; ++i) {
        if (fabs(b[vec->dim - 1]) < tau || isnan(b[vec->dim - 1])) {
            *stt = isnan(b[vec->dim - 1]) ? 14 : 13;

            free_vector(new);
            free(b), free(c), free(d);

            return NULL;
        }
        x[vec->dim - i - 1] =
            (d[vec->dim - i - 1] - c[vec->dim - i - 1] * x[vec->dim - i]) /
            b[vec->dim - i - 1];
    }

    free(b);
    free(c);
    free(d);

    return new;
}

gs_t *new_gs(unsigned dim)
{
    if (!dim)
        return NULL;

    gs_t *new = (gs_t *)malloc(sizeof(gs_t));
    if (!new)
        goto f1;

    new->x = new_vector(dim);
    if (!new->x)
        goto f;

    vector_fill_zero(new->x);
    new->error = NAN;
    new->iter = 0;

    return new;

f:
    free_gs(new);
f1:
    return new;
}

void free_gs(gs_t *gs)
{
    if (!gs)
        return;
    if (gs->x)
        free_vector(gs->x);
    free(gs);
}

void print_gs(gs_t *gs)
{
    if (!gs)
        return;

    printf("Gauss-Seidel solution:\n");
    printf("    ");
    print_vector_shy(gs->x);
    printf("    ");
    printf("Number of iterations: %u\n", gs->iter);
    printf("    ");
    printf("Error of solution: %e\n", gs->error);
}

void lin_gauss_seidel_matrix_tridiag_i(matrix_tridiag_t *mat, vector_t *_b,
                                       vector_t *_x, int *stt)
{
    *stt = -10;
    if (!mat || !_x || !_b)
        return;

    unsigned i = 0;
    unsigned n = mat->dim1;
    double **B = mat->matrix;
    double *x = _x->vector;
    double *b = _b->vector;

    x[0] = (1. / B[0][1]) * (b[0] - B[0][2] * x[1]);
    for (i = 1; i < n - 1; ++i) {
        x[i] =
            (1. / B[i][1]) * (b[i] - B[i][0] * x[i - 1] - B[i][2] * x[i + 1]);
    }
    x[n - 1] = (1. / B[n - 1][1]) * (b[n - 1] - B[n - 1][0] * x[n - 2]);

    *stt = 0;
}

gs_t *solve_gauss_seidel_matrix_tridiag(matrix_tridiag_t *mat, vector_t *vec,
                                        unsigned max_iter, int *stt)
{
    *stt = -10;
    if (!mat || !vec)
        return NULL;

    *stt = 11;
    if (mat->dim1 != vec->dim)
        return NULL;

    *stt = 12;
    if (!matrix_tridiag_is_diagonal_dominant(mat))
        return NULL;

    *stt = -13;
    gs_t *new = new_gs(mat->dim1);
    if (!new)
        return NULL;

    *stt = 0;
    double error = 0.;
    unsigned i = 0;
    for (i = 0; i < max_iter; ++i) {
        lin_gauss_seidel_matrix_tridiag_i(mat, vec, new->x, stt);
        new->iter += 1;
        error = lin_system_error_matrix_tridiag(mat, vec, new->x);
        new->error = error;
        if (error < tau)
            break;
    }

    if (new->iter >= max_iter && new->error > tau)
        *stt = 15;

    new->error = lin_system_error_matrix_tridiag(mat, vec, new->x);

    return new;
}

vector_t *solve_mse_cholesky(matrix_t *mat, vector_t *vec, int *stt)
{
    *stt = -10;
    if (!mat || !vec)
        return NULL;

    *stt = -12;
    matrix_t *trans = matrix_transpose(mat);

    if (!trans)
        return NULL;

    *stt = -12;
    matrix_t *A = matrix_mul_matrix_t(trans);
    if (!A) {
        free_matrix(trans);
        return NULL;
    }

    *stt = -11;
    vector_t *y = matrix_mul_vector(trans, vec);
    if (!y) {
        free_matrix(trans);
        free_matrix(A);
        return NULL;
    }

    vector_t *x = solve_cholesky(A, y, stt);

    free_matrix(trans);
    free_matrix(A);
    free_vector(y);

    return x;
}

ji_t *new_ji(unsigned dim)
{
    if (!dim)
        return NULL;

    ji_t *new = (ji_t *)malloc(sizeof(ji_t));
    if (!new)
        goto f2;

    new->V = new_matrix(dim, dim);
    if (!new->V)
        goto f1;

    new->d = new_vector(dim);
    if (!new->d)
        goto f;

    matrix_fill_identity(new->V);
    vector_fill_zero(new->d);
    new->bmax = NAN;
    new->iter = 0;

    return new;
f:
    free_matrix(new->V);
f1:
    free_ji(new);
f2:
    return new;
}

void free_ji(ji_t *ji)
{
    if (!ji)
        return;
    if (ji->d)
        free_vector(ji->d);
    if (ji->V)
        free_matrix(ji->V);
    free(ji);
}

void print_ji(ji_t *ji)
{
    if (!ji)
        return;

    printf("Eigen-vector and eigen-value decomposition:\n");
    printf("    ");
    printf("Eigen-vectors:\n");
    print_matrix_shy(ji->V);
    printf("    ");
    printf("Eigen-values:\n");
    printf("    ");
    print_vector_shy(ji->d);
    printf("    ");
    printf("Maximmum non-diagonal value on matrix: %e\n", ji->bmax);
    printf("    ");
    printf("Number of iterations: %u\n", ji->iter);
}

ji_t *iterative_jacobi_matrix_eigen(matrix_t *mat, unsigned max_iter, int *stt)
{
    *stt = -10;
    if (!mat)
        return NULL;

    *stt = 11;
    if (mat->dim1 != mat->dim2)
        return NULL;

    *stt = 12;
    if (!matrix_is_sym(mat, tau))
        return NULL;

    *stt = -13;
    ji_t *new = new_ji(mat->dim1);
    if (!new)
        return NULL;

    *stt = -12;
    matrix_t *mat_B = matrix_copy(mat);
    if (!mat_B) {
        free_ji(new);
        return NULL;
    }

    *stt = 0;
    double *bmax = &new->bmax;
    double **B = mat_B->matrix;
    unsigned l = 0, k = 0;
    unsigned n = mat->dim1;

    unsigned i = 0, j = 0;
    double delta = 0.0;
    double t = 0.0;
    double c = 0.0;
    double s = 0.0;
    double bii = 0.0;
    double bjj = 0.0;
    double bij = 0.0;
    for (k = 1; k <= max_iter; ++k) {
        new->iter = k;
        *bmax = matrix_max_diag(mat_B, &i, &j);

        if (*bmax < tau)
            break;

        bii = B[i][i];
        bjj = B[j][j];
        bij = B[i][j];

        delta = (bjj - bii) / (2 * bij);
        t = copysign(1.0, delta) / (fabs(delta) + sqrt(1 + delta * delta));
        c = 1 / sqrt(1 + t * t);
        s = c * t;

        B[i][i] = c * c * bii - 2 * c * s * bij + s * s * bjj;
        B[j][j] = s * s * bii + 2 * c * s * bij + c * c * bjj;
        for (l = 0; l < n; ++l) {
            bii = B[i][l];
            bjj = B[j][l];
            if (l != i && l != j) {
                B[i][l] = c * bii - s * bjj;
                B[j][l] = s * bii + c * bjj;
                B[l][i] = B[i][l];
                B[l][j] = B[j][l];
            }
        }
        B[i][j] = 0.0;
        B[j][i] = 0.0;
        matrix_mul_matrix_givens_i(new->V, i, j, s, c);
    }
    vector_diag_of_matrix_i(new->d, mat_B);

    if (new->iter >= max_iter && new->bmax > tau)
        *stt = 15;

    free_matrix(mat_B);

    return new;
}
