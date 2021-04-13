#ifndef _VECAUX_H
#define _VECAUX_H

#include <math.h>

void __sub_v(double *sink, double *source, unsigned len);
void __sum_v(double *sink, double *source, unsigned len);
void __mul_v(double *sink, double *source, unsigned len);
void __div_v(double *sink, double *source, unsigned len);

void __sub_d(double *vec, unsigned len, double val);
void __sum_d(double *vec, unsigned len, double val);
void __scl_d(double *vec, unsigned len, double val);

void __fnc_v(double *vec, unsigned len, double (*f)(double));

double __dot_product(double *vec1, double *vec, unsigned len);

double __norm_v(double *vec, unsigned len, double p);

double __max_v(double *vec, unsigned len, unsigned *idx);
double __min_v(double *vec, unsigned len, unsigned *idx);

double __max_vabs(double *vec, unsigned len, unsigned *idx);
double __min_vabs(double *vec, unsigned len, unsigned *idx);

#endif
