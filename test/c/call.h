#include "stddef.h"

void callee(int x, float *y) {
    *y = ++x;
}

int caller(void) {
	int x = 0;
	float y = 0;
	double z, w = 9;
	callee(x, &y);
	return y;
}
