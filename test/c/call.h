#include "stddef.h"

int callee(int x, float y) {
    y = ++x;
    return y;
}

int caller(void) {
	int x = 0;
	float y = 0;
	double z, w = 9;
	return callee(x, y);
}
