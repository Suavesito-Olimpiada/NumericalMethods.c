#ifndef STATS_H
#define STATS_H 1

#include <math.h>

#include "vecAux.h"
#include "vector.h"

/* minimum and maximum of a vector */
double min_vector(vector_t *vec);
double max_vector(vector_t *vec);

/* average and variance of a vector */
double avg_vector(vector_t *vec);
double var_vector(vector_t *vec);

#endif /* end of include guard: STATS_H */
