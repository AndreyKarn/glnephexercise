#include <include/libglnsvpos/func.h>
#include <include/libglnsvpos/glnsvpos.h>
#include <include/libglnsvpos/rungekutta.h>
#include <include/libglnsvpos/structures.h>

using namespace std;

uint16_t NT_calc(uint8_t N4, uint16_t T_year, uint16_t year_idx, uint16_t T_month, uint16_t T_day) {
    uint16_t NT;
    uint16_t N42;
    uint16_t Month[12] = {31, 30, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if (T_month < 1 || T_month > 12)
        return 0;

    N42 = (year_idx*31) + N4;
    NT = T_year - 1996 - 4 * (N42 - 1);
    NT *= 365;

    for (uint16_t i = 1; i < T_month; i++)
        NT += Month[i-1];

    return NT + T_day + 1;
}

double GMST_calc(uint8_t N4, uint16_t NT) {
    // Текущая Юлианская дата на 0 часов шкалы МДВ
    double JD0 = 1461 * (N4 - 1) + NT + 2450082.5 - ((NT - 3) / (double)25);
    // Время от эпохи 2000 г 1 января 12 ч (UTC(SU))
    double T_delta = (JD0 - 2451545) / (double)36525;
    // Угол поворота Земли [рад]
    double ERA = 2 * M_PI * ( 0.7790572732640 + 1.00273781191135448 * (JD0 - 2451545));
    // Среднее звездное время по Гринвичу [рад]
    double GMST = ERA;
    GMST += 0.0000000703270726;
    GMST += 0.0223603658710194 * T_delta;
    GMST += 0.0000067465784654 * T_delta * T_delta;
    GMST -= 0.0000000000021332 * T_delta * T_delta * T_delta;
    GMST -= 0.0000000001452308 * T_delta * T_delta * T_delta * T_delta;
    GMST -= 0.0000000000001784 * T_delta * T_delta * T_delta * T_delta * T_delta;

    return GMST;
}

Ephemeris_s CrdTrnsf2Inertial(struct Ephemeris_s Eph, double GMST) {

    struct Ephemeris_s Eph0;

    double Omega_E = 7.2921151467e-5;

    double Theta_Ge = GMST + Omega_E * (Eph.tb - 3 * 60 * 60);

    // Координаты:
    Eph0.X = Eph.X * cos(Theta_Ge) - Eph.Y * sin(Theta_Ge);
    Eph0.Y = Eph.X * sin(Theta_Ge) + Eph.Y * cos(Theta_Ge);
    Eph0.Z = Eph.Z;

    // Скорости:
    Eph0.VX = Eph.VX * cos(Theta_Ge) - Eph.VY * sin(Theta_Ge) - Omega_E * Eph0.Y;
    Eph0.VY = Eph.VX * sin(Theta_Ge) + Eph.VY * cos(Theta_Ge) + Omega_E * Eph0.X;
    Eph0.VZ = Eph.VZ;

    // Ускорения:
    Eph0.AX = Eph.AX * cos(Theta_Ge) - Eph.AY * sin(Theta_Ge);
    Eph0.AY = Eph.AX * sin(Theta_Ge) + Eph.AY * cos(Theta_Ge);
    Eph0.AZ = Eph.AZ;

    // Время
    Eph0.N4 = Eph.N4;
    Eph0.NT = Eph.NT;
    Eph0.tb = Eph.tb;

    return Eph0;
}

void write_struct_Y(struct Y_s *Y_data, uint64_t Size, char *fname) {

    // Запись в файл (для матлаба)
    FILE *file;
    if ((file = fopen(fname,"wb")) == NULL) {
        perror("Error. Problem with file\n");
    }
    else {
        for(uint32_t i = 0; i <= Size; i++) {
            fprintf(file,"%.15e %.15e %.15e %.15e %.15e %.15e\n", Y_data[i].X, Y_data[i].Y, Y_data[i].Z, Y_data[i].VX, Y_data[i].VY, Y_data[i].VZ);
        }
    }
    fclose(file);
}

void read_struct_Y(struct Y_s *Y_data, uint64_t Size, char *fname) {

    // Чтение из файла (из матлаба)

    ifstream file(fname);
    if (file.is_open()) { //Если открытие файла прошло успешно

        string line; //Строчка текста

        uint32_t i;

        for (i = 0; i <= Size; i++) {
            getline(file, line);
            istringstream iss(line);
            iss >> Y_data[i].X;
        }
        for (i = 0; i <= Size; i++) {
            getline(file, line);
            istringstream iss(line);
            iss >> Y_data[i].Y;
        }
        for (i = 0; i <= Size; i++) {
            getline(file, line);
            istringstream iss(line);
            iss >> Y_data[i].Z;
        }
    }
    else perror("Error. Problem with file\n");
}

int add(int a, int b) {
    return mult(a,b) + b;
}
