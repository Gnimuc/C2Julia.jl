#include "stddef.h"

int goto_func(void) {
    int x;
    goto jump;

    jump:
        x += 1;

    return 0;
}
