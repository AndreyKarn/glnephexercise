#ifndef GLNSVPOS_H
#define GLNSVPOS_H

#include <math.h>
#include <iostream>

#include "structures.h"

double GMST_calc(uint8_t N4, uint16_t NT);
Ephemeris_s CrdTrnsf2Inertial(struct Ephemeris_s Eph, double GMST);

int add(int a, int b);


#endif /* #ifndef GLNSVPOS_H */

