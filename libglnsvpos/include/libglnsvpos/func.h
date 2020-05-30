#ifndef FUNC_H
#define FUNC_H

#include <math.h>
#include <iostream>

#include "structures.h"

uint16_t NT_calc(uint8_t N4, uint16_t T_year, uint16_t year_idx, uint16_t T_month, uint16_t T_day);

double GMST_calc(uint8_t N4, uint16_t NT);

Ephemeris_s CrdTrnsf2Inertial(struct Ephemeris_s Eph, double GMST);

void write_struct_Y(struct Y_s *Y_data, uint64_t Size, char *fname);

void read_struct_Y(struct Y_s *Y_data, uint64_t Size, char *fname);

int add(int a, int b);

#endif // FUNC_H
