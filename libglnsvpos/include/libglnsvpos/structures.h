#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iostream>

struct Ephemeris {
    // Time in Gln
    uint8_t   N4;
    uint16_t  NT;
    uint32_t  tb;
    // Coordinates
    double X, Y, Z;
    // Velocity
    double VX, VY,VZ;
    // Acceleration
    double AX, AY, AZ;
};

struct F0 {
    double F1, F2, F3, F4, F5, F6;
};

#endif // STRUCTURES_H
