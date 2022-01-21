#ifndef _CONJ_GRAD_METHOD_INCLUDED_
#define _CONJ_GRAD_METHOD_INCLUDED_

#include "matrix.h"

int check_mat_sym(matrix_t *);

matrix_t * matrix_add(matrix_t *, matrix_t *, int);

matrix_t * scalar_mull(matrix_t *, double);

matrix_t * conj_grad_solver(matrix_t *);

#endif

