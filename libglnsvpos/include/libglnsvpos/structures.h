#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <iostream>

struct Ephemeris {
    // Time in Gln
    uint8_t   N4;
    uint16_t  NT;
    uint32_t  tb;
    // Coordinates
    double X;
    double Y;
    double Z;
    // Velocity
    double VX;
    double VY;
    double VZ;
    // Acceleration
    double AX;
    double AY;
    double AZ;
};

#endif // STRUCTURES_H
