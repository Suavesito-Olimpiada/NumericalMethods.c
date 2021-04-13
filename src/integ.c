#include "integ.h"

double integ_trap(vector_t *y, vector_t *x)
{
    if (!y || !x || y->dim != x->dim)
        return NAN;

    unsigned i = 0;
    double res = 0.0;
    double *f = y->vector, *_x = x->vector;
    for (i = 1; i < x->dim; ++i) {
        res += 0.5 * (f[i] - f[i - 1]) * (_x[i] - _x[i - 1]);
    }

    return res;
}

romberg_t *new_romberg(double (*f)(double), double x0, double xf, unsigned n)
{
    if (!n)
        return NULL;

    romberg_t *new = (romberg_t *)malloc(sizeof(romberg_t));
    if (!new)
        return NULL;

    new->f = f;
    new->x0 = x0;
    new->xf = xf;
    new->n = n;
    new->R = NULL;

    new->R = new_matrix(n + 1, n + 1);
    if (!new->R) {
        free_romberg(new);
        return NULL;
    }
    matrix_fill_zero(new->R);

    return new;
}

void free_romberg(romberg_t *romb)
{
    if (!romb)
        return;
    if (romb->R)
        free_matrix(romb->R);
    free(romb);
}

void print_romberg(romberg_t *romb)
{
    if (!romb)
        return;

    printf("Romberg method solution:\n");
    printf("    ");
    printf("Generated lower matrix:\n");
    print_matrix_shy(romb->R);
    printf("    ");
    printf("a: %.5e\n", romb->x0);
    printf("    ");
    printf("b: %.5e\n", romb->xf);
    printf("    ");
    printf("Partitions: %u\n", ((unsigned)1) << romb->n);
}

static double integ_trap_rec(double (*f)(double), double x0, double xf,
                             double x, unsigned _i)
{
    double h = (xf - x0) / 2.;
    if (_i == 0) {
        return (h * (f(x0) + f(xf)));
    }

    double sum = 0.0;
    h = (xf - x0) / pow(2.0, _i);
    for (unsigned j = 0; j < (unsigned)1 << (_i - 1); ++j) {
        sum += f(x0 + (2. * j - 1.) * h);
    }

    return (.5 * x + h * sum);
}

void solver_romberg(romberg_t *integ)
{
    if (!integ)
        return;

    if (!integ->R)
        return;

    double x0 = integ->x0;
    double xf = integ->xf;
    double (*f)(double) = integ->f;
    double **R = integ->R->matrix;
    unsigned n = integ->n;

    R[0][0] = integ_trap_rec(f, x0, xf, 0.0, 0);
    for (unsigned i = 1; i <= n; ++i) {
        R[i][0] = integ_trap_rec(f, x0, xf, R[i - 1][0], i);
    }

    for (unsigned j = 1; j <= n; ++j) {
        for (unsigned i = j; i <= n; ++i) {
            R[i][j] = R[i][j - 1] + (1. / (pow(4., j) - 1.)) *
                                        (R[i][j - 1] - R[i - 1][j - 1]);
        }
    }
}
