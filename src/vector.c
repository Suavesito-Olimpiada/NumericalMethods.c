#include "vector.h"

vector_t *new_vector(unsigned dim)
{
    vector_t *new = (vector_t *)malloc(sizeof(vector_t));
    if (!new)
        return new;
    new->dim = 0;

    new->vector = (double *)malloc(sizeof(double) * dim);
    if (!(new->vector)) {
        free(new);
        return NULL;
    }

    new->dim = dim;

    return new;
}

void free_vector(vector_t *vec)
{
    if (!vec)
        return;
    free(vec->vector);
    free(vec);
}

vector_t *fread_vector_bin_by_name(char *name)
{
    FILE *pf = fopen(name, "rb");
    if (!pf)
        return NULL;

    unsigned dim = 0;
    fread(&dim, sizeof(unsigned), 1, pf);
    vector_t *new = new_vector(dim);
    if (!new) {
        fclose(pf);
        return NULL;
    }

    fread(new->vector, sizeof(double), dim, pf);
    fclose(pf);

    return new;
}

vector_t *fread_vector_bin(FILE *pf)
{
    if (!pf)
        return NULL;

    unsigned dim = 0;
    fread(&dim, sizeof(unsigned), 1, pf);
    vector_t *new = new_vector(dim);
    if (!new)
        return NULL;

    fread(new->vector, sizeof(double), dim, pf);

    return new;
}

void fwrite_vector_bin_by_name(char *name, vector_t *vec)
{
    if (!vec)
        return;

    FILE *pf = fopen(name, "wb");
    if (!pf)
        return;

    fwrite(&(vec->dim), sizeof(unsigned), 1, pf);
    fwrite(vec->vector, sizeof(double), vec->dim, pf);

    fclose(pf);
}

void fwrite_vector_bin(FILE *pf, vector_t *vec)
{
    if (!vec)
        return;

    if (!pf)
        return;

    fwrite(&(vec->dim), sizeof(unsigned), 1, pf);
    fwrite(vec->vector, sizeof(double), vec->dim, pf);
}

vector_t *fread_vector_text_by_name(char *name)
{
    FILE *pf = fopen(name, "r");
    if (!pf)
        return NULL;

    unsigned dim = 0;
    fscanf(pf, "%u", &dim);
    vector_t *new = new_vector(dim);
    if (!new) {
        fclose(pf);
        return NULL;
    }

    unsigned i = 0;
    double *it = new->vector;
    for (i = 0; i < dim; ++i, ++it) fscanf(pf, "%lf ", it);

    fclose(pf);

    return new;
}

vector_t *fread_vector_text(FILE *pf)
{
    if (!pf)
        return NULL;

    unsigned dim = 0;
    fscanf(pf, "%u ", &dim);
    vector_t *new = new_vector(dim);
    if (!new)
        return NULL;

    unsigned i = 0;
    double *it = new->vector;
    for (i = 0; i < dim; ++i, ++it) fscanf(pf, "%lf", it);

    return new;
}

void fwrite_vector_text_by_name(char *name, vector_t *vec)
{
    if (!vec)
        return;

    FILE *pf = fopen(name, "w");
    if (!pf)
        return;

    fprintf(pf, "%u\n", vec->dim);

    unsigned i = 0;
    double *it = vec->vector;
    for (i = 0; i < vec->dim; ++i, ++it) {
        fprintf(pf, "%lf ", *it);
    }

    fclose(pf);
}

void fwrite_vector_text(FILE *pf, vector_t *vec)
{
    if (!vec)
        return;

    if (!pf)
        return;

    fprintf(pf, "%u\n", vec->dim);

    unsigned i = 0;
    double *it = vec->vector;
    for (i = 0; i < vec->dim; ++i, ++it) {
        fprintf(pf, "%lf ", *it);
    }
}

void print_vector(vector_t *vec)
{
    if (!vec)
        return;

    unsigned i = 0;
    double *it = vec->vector;
    for (i = 0; i < vec->dim; ++i, ++it) printf("%lf ", *it);
    printf("\n");
}

void print_vector_shy(vector_t *vec)
{
    if (!vec)
        return;

    unsigned i = 0;
    double *it = vec->vector;
    if (vec->dim <= 8) {
        for (i = 0; i < vec->dim; ++i, ++it) printf("%1.3e ", *it);
    }
    else {
        for (i = 0; i < 4; ++i) printf("%1.3e ", *(it + i));
        printf("... ");
        for (i = vec->dim - 4; i < vec->dim; ++i) printf("%1.3e ", *(it + i));
    }
    printf("\n");
}

void vector_copy_i(vector_t *dst, vector_t *src)
{
    if (!dst || !src)
        return;

    if (dst->dim != src->dim)
        return;

    unsigned i = 0;
    double *it = src->vector, *it2 = dst->vector;
    for (i = 0; i < src->dim; ++i, ++it, ++it2) *it2 = *it;
}

vector_t *vector_copy(vector_t *vec)
{
    if (!vec)
        return NULL;

    vector_t *new = new_vector(vec->dim);
    if (!new)
        return NULL;

    vector_copy_i(new, vec);

    return new;
}

void vector_fill(vector_t *vec, double val)
{
    if (!vec)
        return;

    unsigned i = 0;
    double *it = vec->vector;
    for (i = 0; i < vec->dim; ++i, ++it) *it = val;
}

void vector_fill_zero(vector_t *vec)
{
    if (!vec)
        return;

    memset(vec->vector, 0, vec->dim * sizeof(double));
}

void vector_fill_one(vector_t *vec)
{
    if (!vec)
        return;

    vector_fill(vec, 1.0);
}

void vector_evaluate_function_over_interval(vector_t *vec, double (*f)(double),
                                            double x_0, double x_n)
{
    if (!vec)
        return;

    unsigned i = 0;
    for (i = 0; i < vec->dim; ++i) {
        vec->vector[i] = f((double)i * ((x_n - x_0) / (vec->dim - 1)) + x_0);
    }
}

void vector_fill_rand(vector_t *vec)
{
    if (!vec)
        return;

    unsigned i = 0;
    double *it = vec->vector;
    for (i = 0; i < vec->dim; ++i, ++it)
        *it = (double)rand() / (double)RAND_MAX;
}

void vector_fill_basis(vector_t *vec, unsigned n)
{
    if (!vec)
        return;

    if (n >= vec->dim)
        return;

    vector_fill_zero(vec);
    vec->vector[n] = 1.0;
}

vector_t *vector_sum_vector(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return NULL;

    if (vec1->dim != vec2->dim)
        return NULL;

    vector_t *new = vector_copy(vec1);
    if (!new)
        return NULL;

    __sum_v(new->vector, vec2->vector, new->dim);
    return new;
}

vector_t *vector_sub_vector(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return NULL;

    if (vec1->dim != vec2->dim)
        return NULL;

    vector_t *new = vector_copy(vec1);
    if (!new)
        return NULL;

    __sub_v(new->vector, vec2->vector, new->dim);
    return new;
}

vector_t *vector_mul_vector(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return NULL;

    if (vec1->dim != vec2->dim)
        return NULL;

    vector_t *new = vector_copy(vec1);
    if (!new)
        return NULL;

    __mul_v(new->vector, vec2->vector, new->dim);
    return new;
}

vector_t *vector_div_vector(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return NULL;

    if (vec1->dim != vec2->dim)
        return NULL;

    vector_t *new = vector_copy(vec1);
    if (!new)
        return NULL;

    __div_v(new->vector, vec2->vector, new->dim);
    return new;
}

vector_t *vector_sum_constant(vector_t *vec, double val)
{
    if (!vec)
        return NULL;

    vector_t *new = vector_copy(vec);
    if (!new)
        return NULL;

    __sum_d(new->vector, new->dim, val);
    return new;
}

vector_t *vector_sub_constant(vector_t *vec, double val)
{
    if (!vec)
        return NULL;

    vector_t *new = vector_copy(vec);
    if (!new)
        return NULL;

    __sub_d(new->vector, new->dim, val);
    return new;
}

vector_t *vector_scl_constant(vector_t *vec, double val)
{
    if (!vec)
        return NULL;

    vector_t *new = vector_copy(vec);
    if (!new)
        return NULL;

    __scl_d(new->vector, new->dim, val);
    return new;
}

vector_t *vector_fnc(vector_t *vec, double (*f)(double))
{
    if (!vec)
        return NULL;

    vector_t *new = vector_copy(vec);
    if (!new)
        return NULL;

    __fnc_v(new->vector, new->dim, f);
    return new;
}

void vector_sum_vector_i(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return;

    if (vec1->dim != vec2->dim)
        return;

    __sum_v(vec1->vector, vec2->vector, vec1->dim);
}

void vector_sub_vector_i(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return;

    if (vec1->dim != vec2->dim)
        return;

    __sub_v(vec1->vector, vec2->vector, vec1->dim);
}

void vector_mul_vector_i(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return;

    if (vec1->dim != vec2->dim)
        return;

    __mul_v(vec1->vector, vec2->vector, vec1->dim);
}

void vector_div_vector_i(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return;

    if (vec1->dim != vec2->dim)
        return;

    __div_v(vec1->vector, vec2->vector, vec1->dim);
}

void vector_sum_constant_i(vector_t *vec, double val)
{
    if (!vec)
        return;

    __sum_d(vec->vector, vec->dim, val);
}

void vector_sub_constant_i(vector_t *vec, double val)
{
    if (!vec)
        return;

    __sub_d(vec->vector, vec->dim, val);
}

void vector_scl_constant_i(vector_t *vec, double val)
{
    if (!vec)
        return;

    __scl_d(vec->vector, vec->dim, val);
}

void vector_fnc_i(vector_t *vec, double (*f)(double))
{
    if (!vec)
        return;

    __fnc_v(vec->vector, vec->dim, f);
}

double vector_dot_product(vector_t *vec1, vector_t *vec2)
{
    if (!vec1 || !vec2)
        return NAN;

    if (vec1->dim != vec2->dim)
        return NAN;

    return __dot_product(vec1->vector, vec2->vector, vec1->dim);
}

double vector_l2_norm(vector_t *vec)
{
    if (!vec)
        return NAN;

    return __norm_v(vec->vector, vec->dim, 2.0);
}

double vector_l1_norm(vector_t *vec)
{
    if (!vec)
        return NAN;

    return __norm_v(vec->vector, vec->dim, 1.0);
}

double vector_linf_norm(vector_t *vec)
{
    if (!vec)
        return NAN;

    return __norm_v(vec->vector, vec->dim, 0.0);
}

double vector_lp_norm(vector_t *vec, double p)
{
    if (!vec)
        return NAN;

    return __norm_v(vec->vector, vec->dim, p);
}
