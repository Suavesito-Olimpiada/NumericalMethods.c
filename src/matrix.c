#include "matrix.h"

matrix_t *new_matrix(unsigned dim1, unsigned dim2)
{
    matrix_t *new = (matrix_t *)malloc(sizeof(matrix_t));
    if (!new)
        return new;
    new->dim1 = 0, new->dim2 = 0;

    new->data = (double *)malloc(sizeof(double) * dim1 * dim2);
    if (!(new->data)) {
        new->matrix = NULL;
        free(new);
        return NULL;
    }

    new->matrix = (double **)malloc(sizeof(double *) * dim1);
    if (!(new->matrix)) {
        free(new->data);
        free(new);
        return NULL;
    }

    unsigned i = 0;
    double **it = new->matrix;
    for (i = 0; i < dim1; ++i, ++it) {
        *it = new->data + i *dim2;
    }

    new->dim1 = dim1, new->dim2 = dim2;

    return new;
}

void free_matrix(matrix_t *mat)
{
    if (!mat)
        return;
    if (!(mat->data))
        return;
    free(mat->data);
    free(mat->matrix);
    free(mat);
}

matrix_t *fread_matrix_bin_by_name(char *name)
{
    FILE *pf = fopen(name, "rb");
    if (!pf)
        return NULL;

    unsigned col = 0, row = 0;
    fread(&row, sizeof(unsigned), 1, pf);
    fread(&col, sizeof(unsigned), 1, pf);
    matrix_t *new = new_matrix(row, col);
    if (!new) {
        fclose(pf);
        return NULL;
    }

    fread(new->data, sizeof(double), col * row, pf);
    fclose(pf);

    return new;
}

matrix_t *fread_matrix_bin(FILE *pf)
{
    if (!pf)
        return NULL;

    unsigned col = 0, row = 0;
    fread(&row, sizeof(unsigned), 1, pf);
    fread(&col, sizeof(unsigned), 1, pf);
    matrix_t *new = new_matrix(row, col);
    if (!new)
        return NULL;

    fread(new->data, sizeof(double), col * row, pf);
    return new;
}

void fwrite_matrix_bin_by_name(char *name, matrix_t *mat)
{
    if (!mat)
        return;

    FILE *pf = fopen(name, "wb");
    if (!pf)
        return;

    fwrite(&(mat->dim1), sizeof(unsigned), 1, pf);
    fwrite(&(mat->dim2), sizeof(unsigned), 1, pf);
    fwrite(mat->data, sizeof(double), mat->dim1 * mat->dim2, pf);

    fclose(pf);
}

void fwrite_matrix_bin(FILE *pf, matrix_t *mat)
{
    if (!mat)
        return;

    if (!pf)
        return;

    fwrite(&(mat->dim1), sizeof(unsigned), 1, pf);
    fwrite(&(mat->dim2), sizeof(unsigned), 1, pf);
    fwrite(mat->data, sizeof(double), mat->dim1 * mat->dim2, pf);
}

matrix_t *fread_matrix_text_by_name(char *name)
{
    FILE *pf = fopen(name, "r");
    if (!pf)
        return NULL;

    unsigned col = 0, row = 0;
    fscanf(pf, "%u %u", &row, &col);
    matrix_t *new = new_matrix(row, col);
    if (!new) {
        fclose(pf);
        return NULL;
    }

    unsigned i = 0, j = 0;
    double *it = new->data;
    for (i = 0; i < row; ++i)
        for (j = 0; j < col; ++j, ++it) fscanf(pf, "%lf", it);

    fclose(pf);

    return new;
}

matrix_t *fread_matrix_text(FILE *pf)
{
    if (!pf)
        return NULL;

    unsigned col = 0, row = 0;
    fscanf(pf, "%u %u", &row, &col);
    matrix_t *new = new_matrix(row, col);
    if (!new)
        return NULL;

    unsigned i = 0, j = 0;
    double *it = new->data;
    for (i = 0; i < row; ++i)
        for (j = 0; j < col; ++j, ++it) fscanf(pf, "%lf", it);

    return new;
}

void fwrite_matrix_text_by_name(char *name, matrix_t *mat)
{
    if (!mat)
        return;

    FILE *pf = fopen(name, "w");
    if (!pf)
        return;

    fprintf(pf, "%u %u\n", mat->dim1, mat->dim2);

    unsigned i = 0, j = 0;
    double *it = mat->data;
    for (i = 0; i < mat->dim1; ++i) {
        for (j = 0; j < mat->dim2; ++j, ++it) fprintf(pf, "%lf ", *it);
        fprintf(pf, "\n");
    }

    fclose(pf);
}

void fwrite_matrix_text(FILE *pf, matrix_t *mat)
{
    if (!mat)
        return;

    if (!pf)
        return;

    fprintf(pf, "%u %u\n", mat->dim1, mat->dim2);

    unsigned i = 0, j = 0;
    double *it = mat->data;
    for (i = 0; i < mat->dim1; ++i) {
        for (j = 0; j < mat->dim2; ++j, ++it) fprintf(pf, "%lf ", *it);
        fprintf(pf, "\n");
    }
}

void print_matrix(matrix_t *mat)
{
    if (!mat)
        return;

    unsigned i = 0, j = 0;
    double *it = mat->data;
    for (i = 0; i < mat->dim1; ++i) {
        for (j = 0; j < mat->dim2; ++j, ++it) printf("%lf ", *it);
        printf("\n");
    }
}

void print_matrix_shy(matrix_t *mat)
{
    if (!mat)
        return;

    unsigned i = 0, j = 0;
    if (mat->dim1 > 8) {
        if (mat->dim2 > 8) {
            for (i = 0; i < 4; ++i) {
                for (j = 0; j < 4; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
                }
                printf(". . .   ");
                for (j = mat->dim2 - 4; j < mat->dim2; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
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
            for (i = mat->dim1 - 4; i < mat->dim1; ++i) {
                for (j = 0; j < 4; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
                }
                printf(". . .   ");
                for (j = mat->dim2 - 4; j < mat->dim2; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
                }
                printf("\n");
            }
        }
        else {
            for (i = 0; i < 5; ++i) {
                for (j = 0; j < mat->dim2; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
                }
                printf("\n");
            }
            for (i = 0; i < 3; ++i) {
                for (j = 0; j < mat->dim2; ++j) printf("%3c.%7c", ' ', ' ');
                printf("\n");
            }
            for (i = mat->dim1 - 5; i < mat->dim1; ++i) {
                for (j = 0; j < mat->dim2; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
                }
                printf("\n");
            }
        }
    }
    else {
        if (mat->dim2 > 10) {
            for (i = 0; i < mat->dim1; ++i) {
                for (j = 0; j < 5; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
                }
                printf("... ");
                for (j = mat->dim2 - 5; j < mat->dim2; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
                }
                printf("\n");
            }
        }
        else {
            for (i = 0; i < mat->dim1; ++i) {
                for (j = 0; j < mat->dim2; ++j) {
                    if (mat->matrix[i][j] < 0.0) {
                        printf("%1.3e ", mat->matrix[i][j]);
                    }
                    else {
                        printf(" %1.3e ", mat->matrix[i][j]);
                    }
                }
                printf("\n");
            }
        }
    }
}

void matrix_copy_i(matrix_t *dst, matrix_t *src)
{
    if (!dst || !src)
        return;

    if (dst->dim1 != src->dim1 || dst->dim2 != src->dim2)
        return;

    unsigned i = 0;
    double *it = src->data, *it2 = dst->data;
    for (i = 0; i < src->dim1 * src->dim2; ++i, ++it, ++it2) *it2 = *it;
}

matrix_t *matrix_copy(matrix_t *mat)
{
    if (!mat)
        return NULL;

    matrix_t *new = new_matrix(mat->dim1, mat->dim2);
    if (!new)
        return NULL;

    matrix_copy_i(new, mat);

    return new;
}

void matrix_fill(matrix_t *mat, double val)
{
    if (!mat)
        return;

    unsigned i = 0;
    double *it = mat->data;
    for (i = 0; i < mat->dim1 * mat->dim2; ++i, ++it) *it = val;
}

void matrix_fill_zero(matrix_t *mat)
{
    if (!mat)
        return;

    memset(mat->data, 0, mat->dim1 * mat->dim2 * sizeof(double));
}

void matrix_fill_rand(matrix_t *mat)
{
    if (!mat)
        return;

    unsigned i = 0;
    double *it = mat->data;
    for (i = 0; i < mat->dim1 * mat->dim2; ++i, ++it)
        *it = (double)rand() / (double)RAND_MAX;
}

void matrix_fill_rand_up(matrix_t *mat)
{
    if (!mat)
        return;

    matrix_fill_zero(mat);
    unsigned i = 0, j = 0;
    for (i = 0; i < mat->dim1; ++i)
        if (i < mat->dim2)
            for (j = i; j < mat->dim2; ++j)
                mat->matrix[i][j] = (double)rand() / (double)RAND_MAX;
}

void matrix_fill_rand_low(matrix_t *mat)
{
    if (!mat)
        return;

    matrix_fill_zero(mat);
    unsigned i = 0, j = 0;
    for (i = 0; i < mat->dim1; ++i)
        if (i < mat->dim2)
            for (j = 0; j <= i; ++j)
                mat->matrix[i][j] = (double)rand() / (double)RAND_MAX;
}

void matrix_fill_identity(matrix_t *mat)
{
    if (!mat)
        return;
    matrix_fill_zero(mat);
    unsigned i = 0;
    for (i = 0; i < ((mat->dim1 < mat->dim2) ? mat->dim1 : mat->dim2); ++i)
        mat->matrix[i][i] = 1.0;
}

void matrix_fill_rand_diag(matrix_t *mat)
{
    if (!mat)
        return;
    matrix_fill_zero(mat);
    unsigned i = 0;
    for (i = 0; i < ((mat->dim1 < mat->dim2) ? mat->dim1 : mat->dim2); ++i)
        mat->matrix[i][i] = (double)rand() / (double)RAND_MAX;
}

void matrix_transpose_i(matrix_t *mat)
{
    if (!mat)
        return;

    if (mat->dim1 != mat->dim2)
        return;

    unsigned i = 0, j = 0;
    double **M = mat->matrix;
    double faux;
    for (j = 0; j < mat->dim2; ++j) {
        for (i = 0; i < j; ++i) {
            faux = M[i][j];
            M[i][j] = M[j][i];
            M[j][i] = faux;
        }
    }
}

matrix_t *matrix_transpose(matrix_t *mat)
{
    if (!mat)
        return NULL;

    matrix_t *new = new_matrix(mat->dim2, mat->dim1);
    if (!new)
        return NULL;

    unsigned i = 0, j = 0;
    for (i = 0; i < mat->dim1; ++i)
        for (j = 0; j < mat->dim2; ++j)
            *(new->data + j * new->dim2 + i) = *(mat->data + i * mat->dim2 + j);

    return new;
}

bool matrix_is_sym(matrix_t *mat, double _tau)
{
    if (!mat)
        return false;

    if (mat->dim1 != mat->dim2)
        return false;

    unsigned i = 0, j = 0;
    double **M = mat->matrix;
    for (j = 0; j < mat->dim2; ++j)
        for (i = 0; i < j; ++i)
            if (fabs(M[i][j] - M[j][i]) > _tau)
                return false;

    return true;
}

matrix_t *matrix_sum_matrix(matrix_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if ((mat1->dim1 != mat2->dim1) || (mat1->dim2 != mat2->dim2))
        return NULL;

    matrix_t *new = matrix_copy(mat1);
    if (!new)
        return NULL;

    __sum_v(new->data, mat2->data, new->dim1 * new->dim2);
    return new;
}

matrix_t *matrix_sub_matrix(matrix_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if ((mat1->dim1 != mat2->dim1) || (mat1->dim2 != mat2->dim2))
        return NULL;

    matrix_t *new = matrix_copy(mat1);
    if (!new)
        return NULL;

    __sub_v(new->data, mat2->data, new->dim1 * new->dim2);
    return new;
}

matrix_t *matrix_div_matrix(matrix_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return NULL;

    if ((mat1->dim1 != mat2->dim1) || (mat1->dim2 != mat2->dim2))
        return NULL;

    matrix_t *new = matrix_copy(mat1);
    if (!new)
        return NULL;

    __div_v(new->data, mat2->data, new->dim1 * new->dim2);
    return new;
}

matrix_t *matrix_sum_constant(matrix_t *mat, double val)
{
    if (!mat)
        return NULL;

    matrix_t *new = matrix_copy(mat);
    if (!new)
        return NULL;

    __sum_d(new->data, new->dim1 * new->dim2, val);
    return new;
}

matrix_t *matrix_sub_constant(matrix_t *mat, double val)
{
    if (!mat)
        return NULL;

    matrix_t *new = matrix_copy(mat);
    if (!new)
        return NULL;

    __sub_d(new->data, new->dim1 * new->dim2, val);
    return new;
}

matrix_t *matrix_scl_constant(matrix_t *mat, double val)
{
    if (!mat)
        return NULL;

    matrix_t *new = matrix_copy(mat);
    if (!new)
        return NULL;

    __scl_d(new->data, new->dim1 * new->dim2, val);
    return new;
}

matrix_t *matrix_fnc(matrix_t *mat, double (*f)(double))
{
    if (!mat)
        return NULL;

    matrix_t *new = matrix_copy(mat);
    if (!new)
        return NULL;

    __fnc_v(new->data, new->dim1 * new->dim2, f);
    return new;
}

void matrix_sum_matrix_i(matrix_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return;

    if ((mat1->dim1 != mat2->dim1) || (mat1->dim2 != mat2->dim2))
        return;

    __sum_v(mat1->data, mat2->data, mat1->dim1 * mat1->dim2);
    return;
}

void matrix_sub_matrix_i(matrix_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return;

    if ((mat1->dim1 != mat2->dim1) || (mat1->dim2 != mat2->dim2))
        return;

    __sub_v(mat1->data, mat2->data, mat1->dim1 * mat1->dim2);
    return;
}

void matrix_div_matrix_i(matrix_t *mat1, matrix_t *mat2)
{
    if (!mat1 || !mat2)
        return;

    if ((mat1->dim1 != mat2->dim1) || (mat1->dim2 != mat2->dim2))
        return;

    __div_v(mat1->data, mat2->data, mat1->dim1 * mat1->dim2);
    return;
}

void matrix_sum_constant_i(matrix_t *mat, double val)
{
    if (!mat)
        return;

    __sum_d(mat->data, mat->dim1 * mat->dim2, val);
    return;
}

void matrix_sub_constant_i(matrix_t *mat, double val)
{
    if (!mat)
        return;

    __sub_d(mat->data, mat->dim1 * mat->dim2, val);
    return;
}

void matrix_scl_constant_i(matrix_t *mat, double val)
{
    if (!mat)
        return;

    __scl_d(mat->data, mat->dim1 * mat->dim2, val);
    return;
}

void matrix_fnc_i(matrix_t *mat, double (*f)(double))
{
    if (!mat)
        return;

    __fnc_v(mat->data, mat->dim1 * mat->dim2, f);
    return;
}

double matrix_nat1_norm(matrix_t *mat)
{
    if (!mat)
        return NAN;

    double max = 0.0, aux = 0.0;
    unsigned i = 0;
    for (i = 0; i < mat->dim1; ++i) {
        aux = __norm_v(mat->data + i * mat->dim2, mat->dim1, 1.0) -
              mat->matrix[i][i];
        if (max < aux)
            max = aux;
    }
    return max;
}

double matrix_natinf_norm(matrix_t *mat)
{
    if (!mat)
        return NAN;

    matrix_t *new = matrix_transpose(mat);
    if (!new)
        return NAN;

    double max = 0.0, aux = 0.0;
    unsigned i = 0;
    for (i = 0; i < new->dim1; ++i) {
        aux = __norm_v(new->data + i * new->dim2, new->dim1, 1.0) -
              mat->matrix[i][i];
        if (max < aux)
            max = aux;
    }

    free_matrix(new);

    return max;
}

double matrix_frobenius_norm(matrix_t *mat)
{
    if (!mat)
        return NAN;

    return __norm_v(mat->data, mat->dim1 * mat->dim2, 2.0);
}

double matrix_max_norm(matrix_t *mat)
{
    if (!mat)
        return NAN;

    unsigned idx;
    return __max_vabs(mat->data, mat->dim1 * mat->dim2, &idx);
}

double matrix_max_diag(matrix_t *mat, unsigned *i, unsigned *j)
{
    if (!mat)
        return NAN;

    if (mat->dim1 < 2 || mat->dim2 < 2)
        return NAN;

    if (!i || !j)
        return NAN;

    double **M = mat->matrix;
    double maxv = fabs(M[0][1]);
    *i = 0;
    *j = 1;
    unsigned _i = 0;
    unsigned _j = 0;
    for (_i = 0; _i < mat->dim1 - 1; ++_i) {
        for (_j = _i + 1; _j < mat->dim2; ++_j) {
            if (maxv < fabs(M[_i][_j])) {
                maxv = fabs(M[_i][_j]);
                *i = _i;
                *j = _j;
            }
        }
    }

    return maxv;
}
