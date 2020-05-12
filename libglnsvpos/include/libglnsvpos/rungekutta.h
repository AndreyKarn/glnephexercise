#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include <math.h>
#include <iostream>

#include "structures.h"

Y_s diffs(double tn , struct Y_s* F);

int RK(struct Ephemeris_s Eph);

int mult(int a, int b);

#endif /* #ifndef RUNGEKUTTA_H */

