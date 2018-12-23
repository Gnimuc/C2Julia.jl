#include "stddef.h"

int loop_stmt(void) {
    int x = 0;

    // for-loop
    // decl inside condition (not C89)
    for ( int i = 0; i < 2; ++i) {
        x += 1;
    }

    // C89 style
    int j;
    for ( j = 0; j < 10; j = j + 1) {
        x = x - 1;
    }

    // while loop
    int k = 0;
    while ( k < 3 ) {
        x = 0;
        k = k + 1;
    }

    // do-while loop
    int m = 3;
    do {
        m = m + 1;
    } while( m < 7 );

    // infinite for-loop
    for( ; ; ) {
        x = 10;
    }

    // infinite for-loop
    for( ; ; ) {
        x = 0;
    }

    return 0;
}
