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
   *xp;

   return 0;
}
