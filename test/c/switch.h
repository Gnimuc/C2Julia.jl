#include "stddef.h"

int switch_func(void) {
    int x = 0;
    char str = 'c';
    switch(str) {
        case 'a' :
            x = 1;
            x = 0;
            break;
            case 'b' :
                if ( x == 2 )
                    x = 7;
        case 'c' :
        default :
            x = 10;
    }
    return 0;
}
