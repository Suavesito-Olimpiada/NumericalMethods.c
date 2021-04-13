#include "polyn.h"

#include <math.h>

polyn_t *new_polyn(unsigned dim) { return new_vector(dim + 1); }

void free_polyn(polyn_t *pol) { free_vector(pol); }

void invert_order_polyn(polyn_t *pol)
{
    if (!pol)
        return;

    unsigned i = 0;
    unsigned n = pol->dim;
    double aux;
    double *p = pol->vector;
    for (i = 0; i < (n >> 1); ++i) {
        aux = p[i];
        p[i] = p[n - i - 1];
        p[n - i - 1] = aux;
    }
}

double evaluate(polyn_t *pol, double x)
{
    if (pol == NULL)
        return NAN;

    unsigned i = 0;
    double res = 0.0;
    for (i = 0; i < pol->dim; ++i) {
        res += pol->vector[i] * pow(x, (double)i);
    }

    return res;
}

double evaluate_horner(polyn_t *pol, double x)
{
    if (pol == NULL)
        return NAN;

    int i = 0;
    double res = 0.0;
    for (i = pol->dim - 1; i >= 0; --i) {
        res = pol->vector[i] + (x * res);
    }

    return res;
}

void print_polyn(polyn_t *pol)
{
    if (!pol)
        return;

    unsigned i = 0;
    for (i = 0; i < pol->dim - 1; ++i) {
        printf("%.3e*x^%d + ", pol->vector[i], i);
    }
    printf("%.3e*x^%d", pol->vector[i], i);
    printf("\n");
}
