#ifndef _CONJ_GRAD_H_
#define _CONJ_GRAD_H_
#include "matrix.h"
int check_mat_sym( matrix_t *mat );
matrix_t *matrix_add( matrix_t *a, matrix_t *b, int alfa );
matrix_t *mnozenie( double a, matrix_t *m );
matrix_t *conj_grad_solver( matrix_t *mat );
#endif
