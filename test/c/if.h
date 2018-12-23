#include "stddef.h"

int if_stmt(void) {
    int a = 0;
    int b = 1;

    // if-then-stmt
    if ( 0 == a ) {
        b = 2;
    }

    // if-else-stmt
    if ( 0 == a ) {
        b = 2;
    } else {
        b = 3;
    }

    // nested if-elseif
    if ( 0 == a ) {
        b = 2;
    } else if ( a < 0 ) {
        b = 3;
    } else {
        b = 4;
    }

    // nested
    if ( 0 == a ) {
        b = 2;
    } else {
        if ( a < 0 ) {
            b = 3;
        } else {
            b = 4;
        }
    }
}
