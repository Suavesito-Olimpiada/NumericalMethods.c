/*  "Integración numérica por el método de Romberg."
 *  autor       : José Joaquín Zubieta Rico
 *  date        : 11-15-2020
 *  subject     : Numerical Methods
 *  description : Integración numérica por el método de Romberg
 * con pivoteo. license     : GPL v3
 * (https://www.gnu.org/licenses/gpl-3.0.en.html)
 */

#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "matrix.h"
#include "vector.h"

double tau;

#include "diffEq.h"

char *argv0;

#include "arg.h"

void usage();

double f(double x, double y);
double y(double x);

int main(int argc, char *argv[])
{
    char *partition_str = NULL;

    ARGBEGIN
    {
        case 'n':
            partition_str = EARGF(usage());
            break;
        default:
            usage();
            break;
    }
    ARGEND;

    unsigned partition = atol(partition_str);

    pvi_t *pvi = new_pvi(f, 4.);
    if(!pvi)
        exit(EXIT_FAILURE);

    double a = 1.;
    double b = 6.;
    rk_t *rk2 = new_rk(pvi, a, b, partition);
    if(!rk2) {
        free_pvi(pvi);
        exit(EXIT_FAILURE);
    }

    vector_t *sol = solve_rk(rk2, 2);
    if(!sol) {
        free_pvi(pvi);
        free_rk(rk2);
        exit(EXIT_FAILURE);
    }

    double *y_ = sol->vector;

    unsigned i = 0;
    double xi = a;
    double error = fabs(y_[0] - y(xi)) / fabs(y(xi));
    double max = error;
    for(i = 0; i < partition; ++i) {
        xi = a + ((double)(i+1) * (b - a)) / (double)(partition);
        error = fabs(y_[i+1] - y(xi)) / fabs(y(xi));
        if (error > max) {
            max = error;
        }
    }

    printf("El error máximo es %.5e.\n", max);

    free_pvi(pvi);
    free_rk(rk2);
    free_vector(sol);

    return 0;
}

double f(double x, double y) { return 4. * x * x - 6. * x + y / x; }

double y(double x) { return -6. * pow(x, 2.) + 8. * x + 2. * pow(x, 3.0); }

void usage()
{
    fprintf(stderr, "usage: %s <-n partition>\n", argv0);
    exit(EXIT_FAILURE);
}
