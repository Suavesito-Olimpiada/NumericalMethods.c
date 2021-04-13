#ifndef _POLYN_H
#define _POLYN_H 1

#include "vector.h"

/* polynomial type */
typedef vector_t polyn_t;

/* polynomial create/free */
polyn_t *new_polyn(unsigned dim);
void free_polyn(polyn_t *pol);

/* invert order of coefficients in polynomial */
void invert_order_polyn(polyn_t *pol);

/* evaluate polynomials */
double evaluate(polyn_t *pol, double x);
double evaluate_horner(polyn_t *pol, double x);

/* input/output fo polynomials */
void print_polyn(polyn_t *pol);

#endif
