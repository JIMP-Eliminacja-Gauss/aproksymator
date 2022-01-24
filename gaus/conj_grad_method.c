
#include <stdlib.h>

#include "conj_grad_method.h"

#ifndef MIN
#define MIN 1e-10
#endif

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
matrix_add(matrix_t *a, matrix_t *b, int alfa) {
    int i = 0;
    int j = 0;
    if (a->rn != b->rn || a->cn != b->cn) {
        return NULL;
    }

    matrix_t *c = make_matrix(a->rn, a->cn);

    if (c == NULL) 
        return NULL;

    for ( i = 0; i < a->rn; i++) {
        for ( j = 0; j < a->cn; j++ ) {
            double a_val = a->e[i * a->cn + j];
            double b_val = b->e[i * b->cn + j];
            put_entry_matrix(c, i, j, a_val + alfa * b_val);
            // int alfa can either be 1 or -1
        }
    }

    return c;
}

matrix_t *
scalar_mull(matrix_t *a, double lambda) {
    int i = 0, j = 0;
    matrix_t *b = copy_matrix(a);

    if (b == NULL) {
        return NULL;
    }

    for (i = 0 ; i < a->rn; i++ ) {
        for (j = 0 ; j < a->cn; j++ ) {
            b->e[ i * a->cn + j ] *= lambda;
        }
    }

    return b;
}


matrix_t * 
conj_grad_solver(matrix_t *mat) {
    int k = 0;
    int j = 0;
    int end = mat->rn;
    double rsold, rsnew;
    double alpha;
    matrix_t *x_mat = make_matrix(mat->rn, 1);
    //testy
    put_entry_matrix(x_mat, 0, 0, 2);
    put_entry_matrix(x_mat, 0, 1, 1);
    //testy
    matrix_t *a_mat = make_matrix(mat->rn, mat->rn);
    matrix_t *b_mat = make_matrix(mat->rn, 1);

    matrix_t *r_mat = NULL;
    matrix_t *p_mat = NULL;
    matrix_t *ap_mat = NULL;
    matrix_t *temp_mat = NULL;
    matrix_t *temp2_mat = NULL;


    if (x_mat == NULL || a_mat == NULL || b_mat == NULL) {
        return NULL;
    }

    for (; k < b_mat->rn; k++) {
        put_entry_matrix(b_mat, k, 0, mat->e[k * mat->cn + mat->cn-1]);
    }


    for (k = 0; k < end; k++) {
        for (j = 0; j < end; j++) {
            put_entry_matrix(a_mat, k, j, mat->e[k * mat->cn + j]);
        }
    }

    // r = b - A * x
    temp_mat = mull_matrix(a_mat, x_mat);
    r_mat = matrix_add(b_mat, temp_mat, -1);
    free_matrix(temp_mat);

    // p = r
    p_mat = copy_matrix(r_mat);

    // rsold = r' * r 
    temp_mat = transpose_matrix(r_mat);
    temp2_mat = mull_matrix(temp_mat, r_mat);
    rsold = temp2_mat->e[0];
    free_matrix(temp_mat);
    free_matrix(temp2_mat);

    for (k = 0; k < end; k++) {
        // Ap = A * p
        if (ap_mat != NULL)
            free_matrix(ap_mat);
        ap_mat = mull_matrix(a_mat, p_mat);

        // alpha = rsold / ( p' * Ap )
        temp_mat = transpose_matrix(p_mat);
        temp2_mat = mull_matrix(temp_mat, ap_mat); 
        alpha = rsold / temp2_mat->e[0]; 
        free_matrix(temp_mat);
        free_matrix(temp2_mat);

        // x = x + alpha * p
        temp_mat = scalar_mull(p_mat, alpha);
        temp2_mat = matrix_add(x_mat, temp_mat, 1);
        free_matrix(x_mat);
        x_mat = copy_matrix(temp2_mat);
        free_matrix(temp_mat);
        free_matrix(temp2_mat);


        // r = r - alpha * Ap
        temp_mat = scalar_mull(ap_mat, alpha);
        temp2_mat = matrix_add(r_mat, temp_mat, -1);
        free_matrix(r_mat);
        r_mat = copy_matrix(temp2_mat);
        free_matrix(temp_mat);
        free_matrix(temp2_mat);

        // rsnew = r' * r
        temp_mat = transpose_matrix(r_mat);
        temp2_mat = mull_matrix(temp_mat, r_mat);
        rsnew = temp2_mat->e[0];
        free_matrix(temp_mat);
        free_matrix(temp2_mat);

        /*if (sqrt(rsnew) < MIN)
            break;*/

        
        // p = r + (rsnew/rsold) * p
        temp_mat = scalar_mull(p_mat, rsnew/rsold);
        free_matrix(p_mat);
        p_mat = matrix_add(r_mat, temp_mat, 1);
        free_matrix(temp_mat);

        // rsold = rsnew    
        rsold = rsnew;
    }
    printf("\n\nPO WSZYSTKIM MACIERZ A\n");
    write_matrix(a_mat, stdout);
    printf("\n\nMACIERZ ROZWIAZAN:\n");
    write_matrix(x_mat, stdout);
    printf("\n\n mat->rn mat->cn %d %d\n", mat->rn, mat->cn );
    j = 0;
    for( k = 0; k < mat->rn*mat->cn; k++ ) {
        if( k % (mat->cn) == 3 ) {
            mat->e[k] = x_mat->e[j];
            j++;
        }
    }
    free_matrix(a_mat);
    free_matrix(b_mat);
    free_matrix(r_mat);
    free_matrix(p_mat);
    free_matrix(ap_mat);

    return x_mat;
}
