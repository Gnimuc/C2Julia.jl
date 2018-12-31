#include "stddef.h"
#include <stdbool.h>

struct foo {
	int x;
};

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

	struct foo f = {x};
	f.x += 1;

	y = z ? 3 : 30;

    return 0;
}
