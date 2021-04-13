#include "stats.h"

double min_vector(vector_t *vec)
{
    if (!vec)
        return NAN;

    unsigned i = 0;
    return __min_v(vec->vector, vec->dim, &i);
}

double max_vector(vector_t *vec)
{
    if (!vec)
        return NAN;

    unsigned i = 0;
    return __max_v(vec->vector, vec->dim, &i);
}
