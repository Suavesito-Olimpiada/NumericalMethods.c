#include "vecAux.h"

void __sub_v(double *sink, double *source, unsigned len)
{
    for (; len--; ++sink, ++source) *sink -= *source;
}

void __sum_v(double *sink, double *source, unsigned len)
{
    for (; len--; ++sink, ++source) *sink += *source;
}

void __mul_v(double *sink, double *source, unsigned len)
{
    for (; len--; ++sink, ++source) *sink *= *source;
}

void __div_v(double *sink, double *source, unsigned len)
{
    for (; len--; ++sink, ++source) *sink /= *source;
}

void __sub_d(double *vec, unsigned len, double val)
{
    for (; len--; ++vec) *vec -= val;
}

void __sum_d(double *vec, unsigned len, double val)
{
    for (; len--; ++vec) *vec += val;
}

void __scl_d(double *vec, unsigned len, double val)
{
    for (; len--; ++vec) *vec *= val;
}

void __fnc_v(double *vec, unsigned len, double (*f)(double))
{
    for (; len--; ++vec) *vec = f(*vec);
}

double __dot_product(double *vec1, double *vec2, unsigned len)
{
    double res = 0.0;
    for (; len--; ++vec1, ++vec2) res += (*vec1) * (*vec2);
    return res;
}

double __norm_v(double *vec, unsigned len, double p)
{
    double res = 0.0;
    if (p == 1.0) {
        for (; len--; ++vec) res += fabs(*vec);
    }
    else if (p == 0.0) {
        res = fabs(*vec);
        for (; len--; ++vec)
            if (fabs(*vec) > res)
                res = *vec;
    }
    else if (p == 2.0) {
        res = __dot_product(vec, vec, len);
        res = sqrt(res);
    }
    else if (p > 0.0) {
        for (; len--; ++vec) res += pow(*vec, p);
        pow(res, 1.0 / p);
    }
    else {
        res = NAN;
    }
    return res;
}

double __max_v(double *vec, unsigned len, unsigned *idx)
{
    double max = *vec;
    for (*idx = 0; *idx < len; ++vec, *idx += 1) {
        if (*vec > max) {
            max = *vec;
        }
    }
    return max;
}

double __min_v(double *vec, unsigned len, unsigned *idx)
{
    unsigned i = *idx = 0;
    double min = *vec;
    for (; len--; ++vec, ++i) {
        if (*vec < min) {
            *idx = i;
            min = *vec;
        }
    }
    return min;
}

double __max_vabs(double *vec, unsigned len, unsigned *idx)
{
    *idx = 0;
    double max = fabs(*vec);
    for (; len--; ++vec, *idx += 1) {
        if (fabs(*vec) > max) {
            max = fabs(*vec);
        }
    }
    return max;
}

double __min_vabs(double *vec, unsigned len, unsigned *idx)
{
    unsigned i = *idx = 0;
    double min = fabs(*vec);
    for (; len--; ++vec, ++i) {
        if (fabs(*vec) < min) {
            *idx = i;
            min = fabs(*vec);
        }
    }
    return min;
}
