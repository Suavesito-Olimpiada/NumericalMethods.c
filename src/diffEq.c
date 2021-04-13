#include "diffEq.h"

pvi_t *new_pvi(double (*f)(double, double), double y0)
{
    if (!f)
        return NULL;

    pvi_t *new = (pvi_t *)malloc(sizeof(pvi_t));
    if (!new)
        return NULL;

    new->y0 = y0;
    new->f = f;

    return new;
}

void free_pvi(pvi_t *pvi)
{
    if (!pvi)
        return;
    free(pvi);
}

void print_pvi(pvi_t *pvi)
{
    if (!pvi)
        return;

    printf("Initial value problem:\n");
    printf("    ");
    printf("Initial condition: %e\n", pvi->y0);
    printf("    Some function...\n");
}

rk_t *new_rk(pvi_t *pvi, double a, double b, unsigned n)
{
    if (!pvi)
        return NULL;

    rk_t *new = (rk_t *)malloc(sizeof(rk_t));
    if (!new)
        return NULL;

    new->a = a;
    new->b = b;
    new->n = n;
    new->pvi = pvi;

    return new;
}

void free_rk(rk_t *rk2)
{
    if (!rk2)
        return;
    free(rk2);
}

void print_rk(rk_t *rk2)
{
    if (!rk2)
        return;

    printf("Runge-Kutta second order method:\n");
    printf("    ");
    printf("Number of iterations: %u\n", rk2->n);
    printf("    ");
    printf("Initial value: %e\n", rk2->a);
    printf("    ");
    printf("Final value: %e\n", rk2->b);
    print_pvi(rk2->pvi);
}

static vector_t *solve_rk2(rk_t *rk2)
{
    if (!rk2)
        return NULL;

    if (!rk2->pvi)
        return NULL;

    unsigned n = rk2->n;

    vector_t *new = new_vector(n + 1);
    if (!new)
        return NULL;

    double a = rk2->a;
    double b = rk2->b;
    double (*f)(double, double) = rk2->pvi->f;
    double *y = new->vector;
    y[0] = rk2->pvi->y0;

    unsigned i = 0;
    double K1 = 0.0, K2 = 0.0;
    double h = (b - a) / (double)(n);
    double xi = 0.0, yi = 0.0;
    for (i = 0; i < n; ++i) {
        yi = y[i];
        xi = a + ((double)(i) * (b - a)) / (double)(n);

        K1 = f(xi, yi);
        K2 = f(xi + h, yi + h * K1);

        y[i + 1] = yi + 0.5 * h * (K1 + K2);
    }

    return new;
}

static vector_t *solve_rk4(rk_t *rk4)
{
    if (!rk4)
        return NULL;

    if (!rk4->pvi)
        return NULL;

    unsigned n = rk4->n;

    vector_t *new = new_vector(n + 1);
    if (!new)
        return NULL;

    double a = rk4->a;
    double b = rk4->b;
    double (*f)(double, double) = rk4->pvi->f;
    double *y = new->vector;
    y[0] = rk4->pvi->y0;

    unsigned i = 0;
    double K1 = 0.0, K2 = 0.0, K3 = 0.0, K4 = 0.0;
    double h = (b - a) / (double)(n);
    double xi = 0.0, yi = 0.0;
    for (i = 0; i < n; ++i) {
        yi = y[i];
        xi = a + ((double)(i) * (b - a)) / (double)(n);

        K1 = h * f(xi, yi);
        K2 = h * f(xi + h / 2., yi + K1 / 2.);
        K3 = h * f(xi + h / 2., yi + K2 / 2.);
        K4 = h * f(xi + h, yi + K3);

        y[i + 1] = yi + (1. / 6.) * (K1 + 2. * K2 + 2. * K3 + K4);
    }

    return new;
}

vector_t *solve_rk(rk_t *rk, unsigned order)
{
    vector_t *res = NULL;
    switch (order) {
        case 2:
            res = solve_rk2(rk);
            break;
        case 4:
            res = solve_rk4(rk);
            break;
        default:
            res = NULL;
            break;
    }
    return res;
}
