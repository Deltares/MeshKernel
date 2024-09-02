#include "rpoly_ak1.h"
#include <stdio.h>

void rpoly(double op [MAX_RPOLY_DEGREE_P1],
           int* Degree,
           double zeror [MAX_RPOLY_DEGREE],
           double zeroi [MAX_RPOLY_DEGREE]) {

  printf ("RPOLY: degree = %i \n", *Degree);

  if (((*Degree) > MAX_RPOLY_DEGREE)||((*Degree) <= 0))
  {
    *Degree = 0;
    return;
  }

  printf ("RPOLY: after degree = %i \n", *Degree);
  printf ("RPOLY: cpeffs = %e %e %e %e %e \n", op[0], op[1], op[2], op[3], op[4]);
  rpoly_ak1 (op, Degree, zeror, zeroi);
  printf ("RPOLY: afterrpoly  degree = %i \n", *Degree);
}
