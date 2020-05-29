#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include <math.h>
#include <iostream>
#include <stdio.h>

#include "structures.h"
#include "func.h"

Y_s* diffs(double tn , struct Y_s Y);

int RK(uint32_t N, double h, struct Y_s* Y) ;

int mult(int a, int b);

#endif /* #ifndef RUNGEKUTTA_H */

