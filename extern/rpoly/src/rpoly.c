#include "rpoly_ak1.h"
#include <stdio.h>

void rpoly(double op [MAX_RPOLY_DEGREE_P1],
           int* Degree,
           double zeror [MAX_RPOLY_DEGREE],
           double zeroi [MAX_RPOLY_DEGREE]) {

  if (((*Degree) > MAX_RPOLY_DEGREE)||((*Degree) <= 0))
  {
    *Degree = 0;
    return;
  }

  rpoly_ak1 (op, Degree, zeror, zeroi);
}
