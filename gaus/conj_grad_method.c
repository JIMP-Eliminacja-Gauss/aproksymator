#include "conj_grad_method.h"

int 
check_mat_sym(matrix_t *mat) {
    if (mat->cn != mat->rn) {
        return -1;
    }

    for (int i = 0; i < mat->rn; i++) {
        for (int j = 0; j < mat->cn; j++) {
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
conj_grad_solver(matrix_t *mat) {
    
}
