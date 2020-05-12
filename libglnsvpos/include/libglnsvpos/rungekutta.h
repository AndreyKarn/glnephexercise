#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include <math.h>
#include <iostream>

#include "structures.h"

F0 diffs(double tn , struct F0* F);

int RK(struct Ephemeris Eph);

int mult(int a, int b);

#endif /* #ifndef RUNGEKUTTA_H */

