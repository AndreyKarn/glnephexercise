#ifndef GLNSVPOS_H
#define GLNSVPOS_H

#include <math.h>
#include <iostream>

#include "structures.h"

double GMST_calc(uint8_t N4, uint16_t NT);
Ephemeris CrdTrnsf2Inertial(struct Ephemeris Eph, double GMST);

int add(int a, int b);


#endif /* #ifndef GLNSVPOS_H */

