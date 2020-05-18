#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iostream>

struct Ephemeris_s {
    // Time in Gln
    uint16_t   N4;
    uint16_t  NT;
    uint32_t  tb;
    // Coordinates
    double X, Y, Z;
    // Velocity
    double VX, VY,VZ;
    // Acceleration
    double AX, AY, AZ;
};

struct Y_s {
    double X, Y, Z, VX, VY, VZ;
};

#endif // STRUCTURES_H
