#include "stddef.h"
#include <stdbool.h>

int operators(void) {
	int x = 0;
    int y = 1;
    bool z = true;

    // unary operators
    +5;
    -5;
    !z;
    (int)z;

    ++y;
    y++;
    --x;
    x--;

    int *xp = &x;
    *xp = 1;

	(bool)(z+1);

	(x == 0) && (y < 1) || (z == 3);

    return 0;
}
