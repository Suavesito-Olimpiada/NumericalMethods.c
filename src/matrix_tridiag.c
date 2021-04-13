#include "matrix_tridiag.h"

matrix_tridiag_t *new_matrix_tridiag(unsigned dim)
{
    return (matrix_tridiag_t *)new_matrix(dim, 3);
}

void free_matrix_tridiag(matrix_tridiag_t *mat)
{
    free_matrix((matrix_t *)mat);
}

matrix_tridiag_t *fread_matrix_tridiag_bin_by_name(char *name)
{
    return (matrix_tridiag_t *)fread_matrix_bin_by_name(name);
}

matrix_tridiag_t *fread_matrix_tridiag_bin(FILE *pf)
{
    return (matrix_tridiag_t *)fread_matrix_bin(pf);
}

void fwrite_matrix_tridiag_bin_by_name(char *name, matrix_tridiag_t *mat)
{
    fwrite_matrix_bin_by_name(name, (matrix_t *)mat);
}

void fwrite_matrix_tridiag_bin(FILE *pf, matrix_tridiag_t *mat)
{
    fwrite_matrix_bin(pf, (matrix_t *)mat);
}

matrix_tridiag_t *fread_matrix_tridiag_text_by_name(char *name)
{
    return (matrix_tridiag_t *)fread_matrix_text_by_name(name);
}

matrix_tridiag_t *fread_matrix_tridiag_text(FILE *pf)
{
    return (matrix_tridiag_t *)fread_matrix_text(pf);
}

void fwrite_matrix_tridiag_text_by_name(char *name, matrix_tridiag_t *mat)
{
    fwrite_matrix_text_by_name(name, (matrix_t *)mat);
}

void fwrite_matrix_tridiag_text(FILE *pf, matrix_tridiag_t *mat)
{
    fwrite_matrix_text(pf, (matrix_t *)mat);
}

void print_matrix_tridiag(matrix_tridiag_t *mat)
{
    if (!mat)
        return;

    unsigned i = 0, j = 0;
    for (j = 1; j < mat->dim1 + 1; ++j) {
        if (j == 1 || j == 2) {
            if (mat->matrix[0][j] < 0.0) {
                printf("%1.3e ", mat->matrix[0][j]);
            }
            else {
                printf(" %1.3e ", mat->matrix[0][j]);
            }
        }
        else
            printf(" %1.3e ", 0.);
    }
    printf("\n");
    for (i = 1; i < mat->dim1 - 1; ++i) {
        for (j = 0; j < mat->dim1; ++j) {
            if (j == i - 1) {
                if (mat->matrix[i][0] < 0.0) {
                    printf("%1.3e ", mat->matrix[i][0]);
                }
                else {
                    printf(" %1.3e ", mat->matrix[i][0]);
                }
            }
            else if (j == i) {
                if (mat->matrix[i][1] < 0.0) {
                    printf("%1.3e ", mat->matrix[i][1]);
                }
                else {
                    printf(" %1.3e ", mat->matrix[i][1]);
                }
            }
            else if (j == i + 1) {
                if (mat->matrix[i][2] < 0.0) {
                    printf("%1.3e ", mat->matrix[i][2]);
                }
                else {
                    printf(" %1.3e ", mat->matrix[i][2]);
                }
            }
            else
                printf(" %1.3e ", 0.);
        }
        printf("\n");
    }
    for (j = 0; j < mat->dim1; ++j) {
        if (j == mat->dim1 - 2) {
            if (mat->matrix[mat->dim1 - 1][0] < 0.) {
                printf("%1.3e ", mat->matrix[mat->dim1 - 1][0]);
            }
            else {
                printf(" %1.3e ", mat->matrix[mat->dim1 - 1][0]);
            }
        }
        else if (j == mat->dim1 - 1) {
            if (mat->matrix[mat->dim1 - 1][1] < 0.) {
                printf("%1.3e ", mat->matrix[mat->dim1 - 1][1]);
            }
            else {
                printf(" %1.3e ", mat->matrix[mat->dim1 - 1][1]);
            }
        }
        else
            printf(" %1.3e ", 0.);
    }
    printf("\n");
}

void print_matrix_tridiag_shy(matrix_tridiag_t *mat)
{
    if (!mat)
        return;

    unsigned i = 0, j = 0;
    if (mat->dim1 > 8) {
        for (j = 1; j < 5; ++j) {
            if (j == 1 || j == 2) {
                if (mat->matrix[0][j] < 0.0) {
                    printf("%1.3e ", mat->matrix[0][j]);
                }
                else {
                    printf(" %1.3e ", mat->matrix[0][j]);
                }
            }
            else
                printf(" %1.3e ", 0.);
        }
        printf(". . .   ");
        for (j = 0; j < 4; ++j) {
            printf(" %1.3e ", 0.);
        }
        printf("\n");
        for (i = 1; i < 4; ++i) {
            for (j = 0; j < 4; ++j) {
                if (j == i - 1) {
                    if (mat->matrix[i][0] < 0.) {
                        printf("%1.3e ", mat->matrix[i][0]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][0]);
                    }
                }
                else if (j == i) {
                    if (mat->matrix[i][1] < 0.) {
                        printf("%1.3e ", mat->matrix[i][1]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][1]);
                    }
                }
                else if (j == i + 1) {
                    if (mat->matrix[i][2] < 0.) {
                        printf("%1.3e ", mat->matrix[i][2]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][2]);
                    }
                }
                else
                    printf(" %1.3e ", 0.);
            }
            printf(". . .   ");
            for (j = 0; j < 4; ++j) {
                printf(" %1.3e ", 0.);
            }
            printf("\n");
        }
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 4; ++j) printf("%3c.%7c", ' ', ' ');
            switch (i) {
                case 0:
                    printf(".       ");
                    break;
                case 1:
                    printf("  .     ");
                    break;
                case 2:
                    printf("    .   ");
                    break;
            }
            for (j = 0; j < 4; ++j) printf("%3c.%7c", ' ', ' ');
            printf("\n");
        }
        printf("\n");
        for (i = mat->dim1 - 4; i < mat->dim1 - 1; ++i) {
            for (j = 0; j < 4; ++j) {
                printf(" %1.3e ", 0.);
            }
            printf(". . .   ");
            for (j = mat->dim1 - 4; j < mat->dim1; ++j) {
                if (j == i - 1) {
                    if (mat->matrix[i][0] < 0.) {
                        printf("%1.3e ", mat->matrix[i][0]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][0]);
                    }
                }
                else if (j == i) {
                    if (mat->matrix[i][1] < 0.) {
                        printf("%1.3e ", mat->matrix[i][1]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][1]);
                    }
                }
                else if (j == i + 1) {
                    if (mat->matrix[i][2] < 0.) {
                        printf("%1.3e ", mat->matrix[i][2]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][2]);
                    }
                }
                else
                    printf(" %1.3e ", 0.);
            }
            printf("\n");
        }
        for (j = 0; j < 4; ++j) {
            printf(" %1.3e ", 0.);
        }
        printf(". . .   ");
        for (j = mat->dim1 - 4; j < mat->dim1; ++j) {
            if (j == mat->dim1 - 2) {
                if (mat->matrix[mat->dim1 - 1][0] < 0.) {
                    printf("%1.3e ", mat->matrix[mat->dim1 - 1][0]);
                }
                else {
                    printf(" %1.3e ", mat->matrix[mat->dim1 - 1][0]);
                }
            }
            else if (j == mat->dim1 - 1) {
                if (mat->matrix[mat->dim1 - 1][1] < 0.) {
                    printf("%1.3e ", mat->matrix[mat->dim1 - 1][1]);
                }
                else {
                    printf(" %1.3e ", mat->matrix[mat->dim1 - 1][1]);
                }
            }
            else
                printf(" %1.3e ", 0.);
        }
        printf("\n");
    }
    else {
        print_matrix_tridiag(mat);
    }
}

matrix_t *matrix_tridiag_to_matrix(matrix_tridiag_t *mat)
{
    if (!mat)
        return NULL;

    matrix_t *new = new_matrix(mat->dim1, mat->dim1);
    if (!new)
        return NULL;

    matrix_fill_zero(new);

    unsigned i = 0;
    for (i = 1; i < mat->dim1 - 1; ++i) {
        new->matrix[i][i - 1] = mat->matrix[i][0];
        new->matrix[i][i] = mat->matrix[i][1];
        new->matrix[i][i + 1] = mat->matrix[i][2];
    }
    new->matrix[0][0] = mat->matrix[0][1];
    new->matrix[0][1] = mat->matrix[0][2];
    new->matrix[new->dim1 - 1][new->dim1 - 2] = mat->matrix[new->dim1 - 1][0];
    new->matrix[new->dim1 - 1][new->dim1 - 1] = mat->matrix[new->dim1 - 1][1];

    return new;
}

matrix_tridiag_t *matrix_to_matrix_tridiag(matrix_t *mat)
{
    if (!mat)
        return NULL;

    matrix_tridiag_t *new = new_matrix_tridiag(mat->dim1);
    if (!new)
        return NULL;

    matrix_tridiag_fill_zero(new);

    unsigned i = 0;
    for (i = 1; i < new->dim1 - 1; ++i) {
        new->matrix[i][0] = mat->matrix[i][i - 1];
        new->matrix[i][1] = mat->matrix[i][i];
        new->matrix[i][2] = mat->matrix[i][i + 1];
    }
    new->matrix[0][1] = mat->matrix[0][0];
    new->matrix[0][2] = mat->matrix[0][1];
    new->matrix[new->dim1 - 1][0] = mat->matrix[new->dim1 - 1][new->dim1 - 2];
    new->matrix[new->dim1 - 1][1] = mat->matrix[new->dim1 - 1][new->dim1 - 1];

    return new;
}

matrix_tridiag_t *matrix_tridiag_copy(matrix_tridiag_t *mat)
{
    return (matrix_tridiag_t *)matrix_copy((matrix_t *)mat);
}

void matrix_tridiag_fill(matrix_tridiag_t *mat, double val)
{
    matrix_fill((matrix_t *)mat, val);
}

void matrix_tridiag_fill_zero(matrix_tridiag_t *mat)
{
    matrix_fill_zero((matrix_t *)mat);
}

void matrix_tridiag_fill_rand(matrix_tridiag_t *mat)
{
    matrix_fill_rand((matrix_t *)mat);
}

void matrix_tridiag_fill_rand_up(matrix_tridiag_t *mat)
{
    if (!mat)
        return;

    unsigned i = 0;
    for (i = 0; i < mat->dim1; ++i) {
        mat->matrix[i][0] = 0.0;
        mat->matrix[i][1] = (double)rand() / (double)RAND_MAX;
        mat->matrix[i][2] = (double)rand() / (double)RAND_MAX;
    }
    mat->matrix[0][0] = 0.0;
    mat->matrix[mat->dim1][2] = 0.0;
}

void matrix_tridiag_fill_rand_low(matrix_tridiag_t *mat)
{
    if (!mat)
        return;

    unsigned i = 0;
    for (i = 0; i < mat->dim1; ++i) {
        mat->matrix[i][0] = (double)rand() / (double)RAND_MAX;
        mat->matrix[i][3] = (double)rand() / (double)RAND_MAX;
        mat->matrix[i][2] = 0.0;
    }
    mat->matrix[0][0] = 0.0;
    mat->matrix[mat->dim1][2] = 0.0;
}

void matrix_tridiag_transpose_i(matrix_tridiag_t *mat)
{
    if (!mat)
        return;

    unsigned i = 0;
    unsigned n = mat->dim1;
    double **A = mat->matrix;
    double aux;
    for (i = 1; i < n; ++i) {
        aux = A[i][0];
        A[i][0] = A[i - 1][2];
        A[i - 1][2] = aux;
    }
}

matrix_tridiag_t *matrix_tridiag_transpose(matrix_tridiag_t *mat)
{
    if (!mat)
        return NULL;

    matrix_tridiag_t *new = matrix_tridiag_copy(mat);
    if (!new)
        return NULL;

    if (new->dim1 == 1)
        return new;

    matrix_tridiag_transpose_i(new);

    return new;
}

bool matrix_tridiag_is_diagonal_dominant(matrix_tridiag_t *mat)
{
    if (!mat)
        return false;

    unsigned i = 0;
    unsigned n = mat->dim1;
    double **A = mat->matrix;
    for (i = 0; i < n; ++i) {
        if (fabs(A[i][1]) < (fabs(A[i][0]) + fabs(A[i][2])))
            return false;
    }

    return true;
}

matrix_tridiag_t *matrix_tridiag_sum_matrix_tridiag(matrix_tridiag_t *mat1,
                                                    matrix_tridiag_t *mat2)
{
    return matrix_sum_matrix((matrix_t *)mat1, (matrix_t *)mat2);
}

matrix_tridiag_t *matrix_tridiag_sub_matrix_tridiag(matrix_tridiag_t *mat1,
                                                    matrix_tridiag_t *mat2)
{
    return matrix_sub_matrix((matrix_t *)mat1, (matrix_t *)mat2);
}

matrix_tridiag_t *matrix_tridiag_div_matrix_tridiag(matrix_tridiag_t *mat1,
                                                    matrix_tridiag_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    matrix_tridiag_t *new = matrix_tridiag_copy(mat1);
    if (!new)
        return NULL;

    unsigned i = 0;
    double *it = new->data + 1, *it2 = mat2->data + 1;
    for (i = 1; i < mat1->dim1; ++i, ++it, ++it2) *it /= *it2;

    return new;
}

matrix_t *matrix_tridiag_sum_matrix(matrix_tridiag_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if (mat1->dim1 != mat2->dim1 || mat1->dim1 != mat2->dim2)
        return NULL;

    matrix_t *new = matrix_tridiag_to_matrix(mat1);
    if (!new)
        return NULL;

    unsigned i = 0, j = 0;
    double *it = mat2->data, *it2 = new->data;
    for (i = 0; i < mat2->dim1; ++i)
        for (j = 0; j < mat2->dim2; ++j, ++it, ++it2) *it2 += *it;

    return new;
}

matrix_t *matrix_tridiag_sub_matrix(matrix_tridiag_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if (mat1->dim1 != mat2->dim1 || mat1->dim1 != mat2->dim2)
        return NULL;

    matrix_t *new = matrix_tridiag_to_matrix(mat1);
    if (!new)
        return NULL;

    unsigned i = 0, j = 0;
    double *it = mat2->data, *it2 = new->data;
    for (i = 0; i < mat2->dim1; ++i)
        for (j = 0; j < mat2->dim2; ++j, ++it, ++it2) *it2 -= *it;

    return new;
}

matrix_tridiag_t *matrix_tridiag_div_matrix(matrix_tridiag_t *mat1,
                                            matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if (mat1->dim1 != mat2->dim1 || mat1->dim1 != mat2->dim2)
        return NULL;

    matrix_tridiag_t *new = matrix_to_matrix_tridiag(mat2);
    if (!new)
        return NULL;

    unsigned i = 0;
    double *it = mat1->data, *it2 = new->data;
    for (i = 1; i < mat1->dim1 * 3 - 1; ++i, ++it, ++it2) {
        *it /= *it2;
    }

    return new;
}

matrix_t *matrix_sub_matrix_tridiag(matrix_t *mat1, matrix_tridiag_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if (mat2->dim1 != mat1->dim1 || mat2->dim1 != mat1->dim2)
        return NULL;

    matrix_t *new = matrix_tridiag_to_matrix(mat2);
    if (!new)
        return NULL;

    unsigned i = 0, j = 0;
    double *it = mat1->data, *it2 = new->data;
    for (i = 0; i < mat1->dim1; ++i) {
        for (j = 0; j < mat1->dim2; ++j, ++it, ++it2) {
            *it2 = *it - *it2;
        }
    }

    return new;
}

matrix_tridiag_t *matrix_tridiag_sum_constant(matrix_tridiag_t *mat, double val)
{
    if (!mat)
        return NULL;

    matrix_tridiag_t *new = matrix_tridiag_copy(mat);
    if (!new)
        return NULL;

    unsigned i = 0;
    double *it = new->data + 1;
    for (i = 1; i < mat->dim1 * 3 - 1; ++i, ++it) *it += val;

    return new;
}

matrix_tridiag_t *matrix_tridiag_scl_constant(matrix_tridiag_t *mat, double val)
{
    if (!mat)
        return NULL;

    matrix_tridiag_t *new = matrix_tridiag_copy(mat);
    if (!new)
        return NULL;

    unsigned i = 0;
    double *it = new->data + 1;
    for (i = 1; i < mat->dim1 * 3 - 1; ++i, ++it) *it *= val;

    return new;
}

double matrix_tridiag_nat1_norm(matrix_tridiag_t *mat)
{
    return matrix_nat1_norm((matrix_t *)mat);
}

double matrix_tridiag_natinf_norm(matrix_tridiag_t *mat)
{
    return matrix_natinf_norm((matrix_t *)mat);
}

double matrix_tridiag_frobenius_norm(matrix_tridiag_t *mat)
{
    return matrix_frobenius_norm((matrix_t *)mat);
}

double matrix_tridiag_max_norm(matrix_tridiag_t *mat)
{
    return matrix_max_norm((matrix_t *)mat);
}
