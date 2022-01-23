#include <stdio.h>
#include <stdlib.h>

#ifndef WYMIAR
#define WYMIAR 3
#endif

int main( int argc, char **argv ) {
    if( argc != 2 ) {
        printf( "zla liczba argumentow, wystarczy jeden argument - nazwa pliku do ktorego wypisac macierz\n" );
        return -1;
    }
    FILE *in = fopen( argv[1], "w" );
    if( in == NULL ) {
        printf( "nie moge pisac do podanego pliku\n" );
        return -2;
    }
    fprintf( in, "%d %d\n", WYMIAR, WYMIAR+1 );
    int macierz[WYMIAR][WYMIAR];
    int macierz_b[WYMIAR];
    int i;
    int j;
    int k=0;
    for( i = 0; i < WYMIAR; i++ )
        macierz_b[i] = /*(*/( (int) rand() )% 100; /*/ RAND_MAX )*100;*/
    for( i = 0; i < WYMIAR; i++ ) {
        for( j = WYMIAR-1; j > k; j-- )
            macierz[i][j] = /*(*/( (int) rand() )%100; /*/ RAND_MAX )*100;*/
        macierz[i][i] = /*(*/( (int) rand() )% 100; /*/ RAND_MAX )*100;*/
        for( j = WYMIAR-1; j > k; j-- )
            macierz[j][i] = macierz[i][j];
        k++;
    }
    for( i = 0; i < WYMIAR; i++ ) {
        for( j = 0; j < WYMIAR; j++ )
            fprintf( in, "%d ", macierz[i][j] );
        fprintf( in, "%d\n", macierz_b[i] );
    }
    fclose(in); 
    return 0;
}

