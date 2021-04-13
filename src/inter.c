#include "inter.h"

polyn_t *approx_polyn(vector_t *_x, vector_t *_y, unsigned n, int *stt)
{
    *stt = -10;
    if (!_x || !_y)
        return NULL;

    *stt = 11;
    if (_x->dim != _y->dim)
        return NULL;

    *stt = -12;
    matrix_t *new = new_matrix(_x->dim, n + 1);
    if (!new)
        return NULL;

    unsigned i = 0, j = 0;
    double *x = _x->vector;
    double **A = new->matrix;
    for (i = 0; i < new->dim1; ++i) {
        for (j = 0; j < new->dim2 - 1; ++j) {
            A[i][j] = pow(x[i], (double)(n - j));
        }
        A[i][new->dim2 - 1] = 1.0;
    }

    vector_t *pol = solve_mse_cholesky(new, _y, stt);

    free_matrix(new);

    return (polyn_t *)pol;
}

matrix_t *dif_div(matrix_t *pts)
{
    if (pts == NULL || pts->dim2 != 2)
        return NULL;

    matrix_t *new = new_matrix(pts->dim1, pts->dim1);
    if (!new)
        return NULL;

    matrix_fill_zero(new);

    unsigned i = 0, j = 0;
    for (i = 0; i < pts->dim1; ++i) new->matrix[i][0] = pts->matrix[i][1];

    double **A = new->matrix, **x = pts->matrix;
    for (j = 1; j < pts->dim1; ++j)
        for (i = 0; i < pts->dim1 - j; ++i)
            A[i][j] = (A[i + 1][j - 1] - A[i][j - 1]) / (x[i + j][0] - x[i][0]);

    return new;
}

double evaluate_matrix_inter_dif_div(matrix_t *coef, matrix_t *val, double x)
{
    unsigned i = 0, j = 0;
    double res = 0.0, aux = 1.0;
    for (i = 0; i < coef->dim1; ++i) {
        aux = 1.0;
        for (j = 0; j < i; ++j) {
            aux *= (x - val->matrix[j][0]);
        }
        res += coef->matrix[0][i] * aux;
    }

    return res;
}

splin3_t *new_splin3(matrix_t *xf)
{
    if (!xf)
        return NULL;

    if (xf->dim2 != 2)
        return NULL;

    splin3_t *new = (splin3_t *)malloc(sizeof(splin3_t));
    if (!new)
        return NULL;

    new->n = xf->dim1 - 1;
    new->x = NULL;
    new->f = NULL;
    new->M = NULL;

    new->x = vector_col_of_matrix(xf, 0);
    if (!new->x) {
        free_splin3(new);
        return NULL;
    }

    new->f = vector_col_of_matrix(xf, 1);
    if (!new->f) {
        free_splin3(new);
        return NULL;
    }

    return new;
}

void free_splin3(splin3_t *splin)
{
    if (!splin)
        return;
    if (splin->x)
        free_vector(splin->x);
    if (splin->f)
        free_vector(splin->f);
    if (splin->M)
        free_vector(splin->M);
    free(splin);
}

void solve_splin3(splin3_t *splin, double d0, double dn, int *stt)
{
    *stt = -10;
    if (!splin)
        return;

    *stt = -12;
    matrix_tridiag_t *mat = new_matrix_tridiag(splin->n + 1);
    if (!mat)
        return;
    matrix_fill_zero(mat);

    *stt = -11;
    vector_t *vec = new_vector(splin->n + 1);
    if (!vec) {
        free_matrix_tridiag(mat);
        return;
    }
    vector_fill_zero(vec);

    double hi = 0.0, hi1 = 0.0;
    double mu = 0.0, lambda = 0.0;
    double *x = splin->x->vector;
    double *f = splin->f->vector;
    double *d = vec->vector;
    double **m = mat->matrix;
    unsigned n = splin->n;
    unsigned i = 0;

    m[0][1] = 2.0;
    m[0][2] = 1.0;
    d[0] = d0;
    for (i = 1; i < n; ++i) {
        hi = (x[i] - x[i - 1]);
        hi1 = (x[i + 1] - x[i]);
        lambda = (hi1) / (hi + hi1);
        mu = 1.0 - lambda;

        m[i][0] = mu;
        m[i][1] = 2.0;
        m[i][2] = lambda;
        d[i] = (6.0 / (hi + hi1)) *
               ((f[i + 1] - f[i]) / hi1 - (f[i] - f[i - 1]) / hi);
    }
    m[n][0] = 1.0;
    m[n][1] = 2.0;
    d[n] = dn;

    splin->M = solve_matrix_tridiag(mat, vec, stt);

    free_matrix_tridiag(mat);
    free_vector(vec);

    return;
}

void solve_nat_splin3(splin3_t *splin, int *stt)
{
    *stt = -10;
    if (!splin)
        return;

    solve_splin3(splin, 0, 0, stt);
}

double eval_splin3(splin3_t *splin, double _x)
{
    if (!splin)
        return NAN;

    if (!splin->M)
        return NAN;

    double *x = splin->x->vector;
    double *f = splin->f->vector;
    double *M = splin->M->vector;

    unsigned n = splin->n;
    if (_x < x[0] || _x > x[n])
        return NAN;

    unsigned i = 0;
    for (i = 0; i < n; ++i) {
        if (x[i] < _x && _x < x[i + 1])
            break;
    }
    ++i;

    double h = x[i] - x[i - 1];
    double c = f[i - 1] - M[i - 1] * (h * h / 6.0);
    double C = (f[i] - f[i - 1]) / h - (h / 6.0) * (M[i] - M[i - 1]);

    double Sx = (M[i - 1] / (6.0 * h)) * pow((x[i] - _x), 3.0) +
                (M[i] / (6.0 * h)) * pow((_x - x[i - 1]), 3.0) +
                C * (_x - x[i - 1]) + c;

    return Sx;
}
