#include "conj_grad_method.h"
#include "matrix.h"

#ifndef MIN
#define MIN 1e-9
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

    for ( ; i < a->rn; i++) {
        for ( ; j < a->cn; j++ ) {
            double a_val = a->e[i * a->cn + j];
            double b_val = b->e[i * b->cn + j];
            add_to_entry_matrix(c, i, j, a_val + alfa * b_val);
            // int alfa can either be 1 or -1
        }
    }

    return c;
}


matrix_t * 
conj_grad_solver(matrix_t *mat) {
    int k = 0;
    int j = 0;
    int end = mat->rn;
    double rsold;
    double alfa;
    double beta;
    matrix_t *x_mat = make_matrix(mat->cn, 1);
    matrix_t *r_mat = make_matrix(mat->cn, 1);
    matrix_t *p_mat = NULL;
    matrix_t *a_mat = make_matrix(mat->rn, mat->rn);
    matrix_t *b_mat = make_matrix(mat->rn, 1);
    matrix_t *tmp_r = NULL;
    matrix_t *tmp_p = NULL;

    if (x_mat == NULL || r_mat == NULL || a_mat == NULL || b_mat == NULL) {
        return NULL;
    }

    for (; k < b_mat->rn; k++) {
        put_entry_matrix(b_mat, k, 0, mat->e[k * mat->cn + mat->cn]);
    }

    for (k = 0; k < end; k++) {
        for (j = 0; j < end; j++) {
            put_entry_matrix(a_mat, k, j, mat->e[k * mat->cn + j]);
        }
    }

    r_mat = matrix_add(b_mat, mull_matrix(a_mat, x_mat), -1);
    p_mat = copy_matrix(r_mat);
    rsold = (mull_matrix(transpose_matrix(r_mat), r_mat))->e[0];
    alfa = rsold / mull_matrix( mull_matrix( transpose_matrix(p_mat), a_mat ), p_mat )->e[0];
    x_mat = mnozenie( alfa, p_mat );
    tmp_r = r_mat;
    tmp_p = p_mat;


    for (k = 0; k < end; k++) {
        r_mat = odejmowanie( tmp_r, mull_matrix( mnozenie( alfa, a_mat ), p_mat ) );
        //x_mat = x_mat + alfa*p1;
        beta = mull_matrix( transpose_matrix(r_mat), r_mat )->e[0] / mull_matrix( transpose_matrix(tmp_r), tmp_r )->e[0];
        p_mat = dodawanie( r_mat, mnozenie( beta, tmp_p ) );
        free(tmp_p);
        free(tmp_r);
        tmp_p = p_mat;
        tmp_r = r_mat;
        alfa = mull_matrix( transpose_matrix(r_mat), r_mat )->e[0] / mull_matrix( mull_matrix( transpose_matrix(p_mat), a_mat ), p_mat )->e[0];
        x_mat = dodawanie( x_mat, mnozenie( alfa, p_mat ) );
    }

    return NULL;
}
