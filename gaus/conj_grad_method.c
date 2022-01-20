#include "conj_grad_method.h"

int 
check_mat_sym(matrix_t *mat) {
    int i, j;

    if (mat->cn != mat->rn) {
        return -1;
    }

    for (i = 0; i < mat->rn; i++) {
        for (j = 0; j < mat->cn; j++) {
            if (mat->e[i * mat->cn + j] != mat->e[j * mat->cn + j]) {
                return 0;
                // meaning matrix is not symmetrical
            }
        }
    }

    return 1;
    // matrix is indeed symmetrical
}

matrix_t *
matrix_substract(matrix_t *a, matrix_t *b) {
    return NULL;

}


matrix_t * 
conj_grad_solver(matrix_t *mat) {
    int k = 0;
    int end = mat->rn;
    matrix_t *x_mat = make_matrix(mat->cn, 1);
    matrix_t *r_mat = make_matrix(mat->cn, 1);
    matrix_t *p_mat = copy_matrix(r_mat);


    for (; k < end; k++) {
        ;
    }

    return NULL;
}
