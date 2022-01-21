#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matrix.h"
#include "conj_grad_method.h"
int main( int argc, char **argv ) {
    FILE *in;
    if( argc != 2 )
        return -1;
    in = fopen( argv[1], "r" );
    matrix_t *m = read_matrix(in);    
    printf("\nMacierz:\n");
    write_matrix( m, stdout );
    printf("\nPo rozwiazaniu:\n");
    write_matrix( conj_grad_solver(m), stdout );
    return 0;
}
