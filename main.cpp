#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/timeb.h>
#include <time.h>

using namespace std;

// Исходные данные
#define CoordinateX 23036950.68
#define CoordinateY -9091173.34
#define CoordinateZ 6041059.08

#define VelocityX 755.86033
#define VelocityY -358.52718
#define VelocityZ -3447.90649

#define AccelerationX 0.0000056
#define AccelerationY 0.0000000
#define AccelerationZ -0.0000028

#define Time_Year 2020
#define Time_Month 2
#define Time_Day 10
#define Time_Hour 13
#define Time_Minutes 45
#define Time_Seconds 18


struct Ephemeris_type {
  uint16_t  N4;     // Time in Gln
  uint16_t  NT;
  uint32_t  tb;
  double X, Y, Z,    // Coordinates
         VX, VY,VZ, // Velocity
         AX, AY, AZ; // Acceleration
};

struct Y_type {
  double X, Y, Z, VX, VY, VZ;
};


//#define _USE_MATH_DEFINES
#include <math.h> // число pi
#ifndef M_PI
 #define M_PI 3.14159265358979323846
#endif
using namespace std;

uint16_t NT_Calculate(uint8_t N4, uint16_t Tyear, uint16_t year_idx, uint16_t Tmonth, uint16_t Tday) {
    uint16_t NT, N42;
    uint8_t Mon[12] = {31, 30, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if (Tmonth < 1 || Tmonth > 12) {
       return 0;
    }

    N42 = (year_idx*31) + N4;
    NT = Tyear - 1996 - 4 * (N42 - 1);
    NT *= 365;

    for (uint16_t i = 1; i < Tmonth; i++) {
      NT += Mon[i-1];
    }

    return (NT + Tday + 1);
}

double GMST_Calculate(uint8_t N4, uint16_t NT) {
    // Текущая Юлианская дата на 0 часов шкалы МДВ
    double JD0 = 1461 * (N4 - 1) + NT + 2450082.5 - ((NT - 3) / 25.0);
    // Время от эпохи 2000 г 1 января 12 ч (UTC(SU))
    double Td = (JD0 - 2451545) / 36525.0;
    // Угол поворота Земли [рад]
    double ERA = 2 * M_PI * ( 0.7790572732640 + 1.00273781191135448 * (JD0 - 2451545));
    // Среднее звездное время по Гринвичу [рад]
    double GMST = ERA;
    GMST += 0.0000000703270726;
    GMST += 0.0223603658710194 * Td;
    GMST += 0.0000067465784654 * Td * Td;
    GMST -= 0.0000000000021332 * Td * Td * Td;
    GMST -= 0.0000000001452308 * Td * Td * Td * Td;
    GMST -= 0.0000000000001784 * Td * Td * Td * Td * Td;

   return GMST;
}


void CrdTrnsfToInertial(struct Ephemeris_type Ephemeris, double GMST, struct Ephemeris_type &Ephemeris0) {

    double OmegaE = 7.2921151467e-5,
           ThetaGe = GMST + OmegaE * (Ephemeris.tb - 3*60*60);

    // Координаты:
    Ephemeris0.X = Ephemeris.X * cos(ThetaGe) - Ephemeris.Y * sin(ThetaGe);
    Ephemeris0.Y = Ephemeris.X * sin(ThetaGe) + Ephemeris.Y * cos(ThetaGe);
    Ephemeris0.Z = Ephemeris.Z;

    // Скорости:
    Ephemeris0.VX = Ephemeris.VX * cos(ThetaGe) - Ephemeris.VY * sin(ThetaGe) - OmegaE * Ephemeris0.Y;
    Ephemeris0.VY = Ephemeris.VX * sin(ThetaGe) + Ephemeris.VY * cos(ThetaGe) + OmegaE * Ephemeris0.X;
    Ephemeris0.VZ = Ephemeris.VZ;

    // Ускорения:
    Ephemeris0.AX = Ephemeris.AX * cos(ThetaGe) - Ephemeris.AY * sin(ThetaGe);
    Ephemeris0.AY = Ephemeris.AX * sin(ThetaGe) + Ephemeris.AY * cos(ThetaGe);
    Ephemeris0.AZ = Ephemeris.AZ;

    // Время
    Ephemeris0.N4 = Ephemeris.N4;
    Ephemeris0.NT = Ephemeris.NT;
    Ephemeris0.tb = Ephemeris.tb;
}

void SaveData(struct Y_type *Y_data, uint64_t Size, char *fname) {
  // Сохранение данных (для матлаба)
   FILE *fs;
   if ((fs = fopen(fname,"wb")) == NULL) {
     printf("Error. File: %s, Line: %d\n", __FILE__, __LINE__);
   } else {
      for(uint32_t i = 0; i <= Size; i++) {
        fprintf(fs, "%.15e %.15e %.15e %.15e %.15e %.15e\n", Y_data[i].X, Y_data[i].Y, Y_data[i].Z, Y_data[i].VX, Y_data[i].VY, Y_data[i].VZ);
      }
      fclose(fs);
   }
}

Y_type* DiffsY(struct Y_type Y)
{
    double J02 = 1082625.75e-9, // зональный гармонический коэффициент второй степени, характеризующий полярное сжатие Земли
     GM = 398600441.8e6, // геоцентрическая константа гравитационного поля Земли с учетом атмосферы, [м3/c2]
     a_e = 6378136; // большая полуось общеземного эллипсоида, [м]

    double crdX = Y.X, crdY = Y.Y, crdZ = Y.Z;
    double r = sqrt(crdX * crdX + crdY * crdY + crdZ * crdZ);

    double GM0 = GM / (r * r);
    double Rho = a_e / r;
    double crdX0 = crdX / r, crdY0 = crdY / r, crdZ0 = crdZ / r;

    struct Y_type *dY = new struct Y_type;
    // Дифуры
    dY->X = Y.VX;
    dY->Y = Y.VY;
    dY->Z = Y.VZ;

    dY->VX = - GM0 * crdX0 - 1.5 * J02 * GM0 * crdX0 * Rho * Rho * (1.0 - 5.0 * crdZ0 * crdZ0);
    dY->VY = - GM0 * crdY0 - 1.5 * J02 * GM0 * crdY0 * Rho * Rho * (1.0 - 5.0 * crdZ0 * crdZ0);
    dY->VZ = - GM0 * crdZ0 - 1.5 * J02 * GM0 * crdZ0 * Rho * Rho * (3.0 - 5.0 * crdZ0 * crdZ0);

    return dY;
}

int RK(uint32_t N, double h, struct Y_type* Y) {

    if (N == 0 || fabs(h)< 0.0000001 ) return 0;

    struct Y_type *k1, *k2, *k3, *k4, *knextstep;
    struct Y_type Y2, Y3, Y4;

    for (uint32_t k = 1; k < N; k++)  {
        k1 = new struct Y_type;
        k2 = new struct Y_type;
        k3 = new struct Y_type;
        k4 = new struct Y_type;
        knextstep = new struct Y_type;

        k1 = DiffsY(Y[k-1]);

        Y2.X = Y[k-1].X + h * k1->X / 2.0;
        Y2.Y = Y[k-1].Y + h * k1->Y / 2.0;
        Y2.Z = Y[k-1].Z + h * k1->Z / 2.0;
        Y2.VX = Y[k-1].VX + h * k1->VX / 2.0;
        Y2.VY = Y[k-1].VY + h * k1->VY / 2.0;
        Y2.VZ = Y[k-1].VZ + h * k1->VZ / 2.0;

        k2 = DiffsY(Y2);

        Y3.X = Y[k-1].X + h * k2->X / 2.0;
        Y3.Y = Y[k-1].Y + h * k2->Y / 2.0;
        Y3.Z = Y[k-1].Z + h * k2->Z / 2.0;
        Y3.VX = Y[k-1].VX + h * k2->VX / 2.0;
        Y3.VY = Y[k-1].VY + h * k2->VY / 2.0;
        Y3.VZ = Y[k-1].VZ + h * k2->VZ / 2.0;

        k3 = DiffsY (Y3);

        Y4.X = Y[k-1].X + h * k3->X;
        Y4.Y = Y[k-1].Y + h * k3->Y;
        Y4.Z = Y[k-1].Z + h * k3->Z;
        Y4.VX = Y[k-1].VX + h * k3->VX;
        Y4.VY = Y[k-1].VY + h * k3->VY;
        Y4.VZ = Y[k-1].VZ + h * k3->VZ;

        k4 = DiffsY (Y4);

        knextstep->X = h / 6.0 * ( k1->X + 2.0 * k2->X + 2.0 * k3->X + k4->X );
        knextstep->Y = h / 6.0 * ( k1->Y + 2.0 * k2->Y + 2.0 * k3->Y + k4->Y );
        knextstep->Z = h / 6.0 * ( k1->Z + 2.0 * k2->Z + 2.0 * k3->Z + k4->Z );
        knextstep->VX = h / 6.0 * ( k1->VX + 2.0 * k2->VX + 2.0 * k3->VX + k4->VX );
        knextstep->VY = h / 6.0 * ( k1->VY + 2.0 * k2->VY + 2.0 * k3->VY + k4->VY );
        knextstep->VZ = h / 6.0 * ( k1->VZ + 2.0 * k2->VZ + 2.0 * k3->VZ + k4->VZ );

        Y[k].X = Y[k-1].X + knextstep->X;
        Y[k].Y = Y[k-1].Y + knextstep->Y;
        Y[k].Z = Y[k-1].Z + knextstep->Z;
        Y[k].VX = Y[k-1].VX + knextstep->VX;
        Y[k].VY = Y[k-1].VY + knextstep->VY;
        Y[k].VZ = Y[k-1].VZ + knextstep->VZ;

        delete k1;
        delete k2;
        delete k3;
        delete k4;
        delete knextstep;
    }
    return 0;
}


uint64_t GlonassVPos(bool RK_valid, double h) {

    // Class Ephemeris
    struct Ephemeris_type Ephemeris;
    // Coordinates
    Ephemeris.X = CoordinateX;
    Ephemeris.Y = CoordinateY;
    Ephemeris.Z = CoordinateZ;

    // Velocity
    Ephemeris.VX = VelocityX;
    Ephemeris.VY = VelocityY;
    Ephemeris.VZ = VelocityZ;

    // Acceleration
    Ephemeris.AX = AccelerationX;
    Ephemeris.AY = AccelerationY;
    Ephemeris.AZ = AccelerationZ;

    uint16_t Time_year = Time_Year;
    uint8_t  Time_month = Time_Month,
             Time_day = Time_Day,
             Time_hour = Time_Hour,
             Time_minutes = Time_Minutes,
             Time_seconds = Time_Seconds;

    uint16_t year_idx = 0;

    // Time in Gln
    Ephemeris.N4 = ((Time_year - 1996) / 4) + 1;
    while (Ephemeris.N4 > 31) { // Учет 5-битности N4
        Ephemeris.N4 -= 31;
        year_idx++;
    }

    Ephemeris.NT = NT_Calculate( Ephemeris.N4, Time_year, year_idx, Time_month, Time_day );

    Ephemeris.tb = Time_seconds + Time_minutes*60 + Time_hour*60*60 + 10800;
    if (Ephemeris.tb >= 24*60*60) {
        Ephemeris.tb -= 24*60*60;
        Ephemeris.NT++;
        if (Ephemeris.NT >= 1462) Ephemeris.N4++;
    }

    double GMST = GMST_Calculate( Ephemeris.N4, Ephemeris.NT );

    struct Ephemeris_type Ephemeris0;
    CrdTrnsfToInertial(Ephemeris, GMST, Ephemeris0);

    uint32_t tn = Ephemeris.tb; // Текущее время
    uint32_t Toe = (12+3)*60*60; // Начальное время
    uint32_t Tof = (24+3)*60*60; // Конечное время

    uint64_t N2inc = (Tof - tn) / h; // Количество отcчетов для времени большего текущего Ephemeris.tb
    uint64_t N2dec = (tn - Toe) / h; // Количество отcчетов для времени меньшего текущего Ephemeris.tb
    uint64_t N = N2inc + N2dec; // Общее число отсчетов

    // Численное интегрирования для времени меньшего текущего Ephemeris.tb
    struct Y_type *Ydec;
    Ydec = new struct Y_type[N2dec];
    Ydec[0].X = Ephemeris0.X;
    Ydec[0].Y = Ephemeris0.Y;
    Ydec[0].Z = Ephemeris0.Z;
    Ydec[0].VX = Ephemeris0.VX;
    Ydec[0].VY = Ephemeris0.VY;
    Ydec[0].VZ = Ephemeris0.VZ;

    if (RK_valid) RK( N2dec, -h, Ydec);

    // Численное интегрирования для времени большего текущего Ephemeris.tb
    struct Y_type *Yinc;
    Yinc = new struct Y_type[N2inc];
    Yinc[0].X = Ephemeris0.X;
    Yinc[0].Y = Ephemeris0.Y;
    Yinc[0].Z = Ephemeris0.Z;
    Yinc[0].VX = Ephemeris0.VX;
    Yinc[0].VY = Ephemeris0.VY;
    Yinc[0].VZ = Ephemeris0.VZ;

    if (RK_valid) RK( N2inc, h, Yinc);

    // Формирование выходных данных
    struct Y_type *Yout;
    Yout = new struct Y_type[N];
    uint32_t i;
    for (i = 0; i < N2dec; i++) {
        Yout[i] = Ydec[N2dec-i-1];
    }
    for (i = 0; i < N2inc; i++) {
        Yout[i+N2dec] = Yinc[i];
    }

    // Учет ускорений
    double tau = (double)Toe - (double)tn;
    for (i = 0; i < N; i++) {
      double tau2 =tau * tau;
      Yout[i].X += Ephemeris0.AX * tau2 / 2.0;
      Yout[i].Y += Ephemeris0.AY * tau2 / 2.0;
      Yout[i].Z += Ephemeris0.AZ * tau2 / 2.0;

      Yout[i].VX += Ephemeris0.AX * tau;
      Yout[i].VY += Ephemeris0.AY * tau;
      Yout[i].VZ += Ephemeris0.AZ * tau;

      tau += h;
    }

    SaveData(Yout, N, (char*)"glonass_data.txt");

    // Очищение памяти
    delete []Yout;
    delete []Ydec;
    delete []Yinc;

    return N;
}

int main() {
/*  бла бла   */
    setlocale(LC_ALL, "Russian");

    cout << "Начало подсчета" << endl;

    struct timeb tb;
    ::ftime(&tb);
    uint64_t st=1000ull*(uint64_t)tb.time+(uint64_t)tb.millitm;
    string msg;
    uint64_t rc=GlonassVPos(0, 0.1);

    ::ftime(&tb);
    uint64_t et=1000ull*(uint64_t)tb.time+(uint64_t)tb.millitm;
    cout << (rc>0?"Подсчет закончен успешно":"Подсчет закончен с ошибкой ") << " " << (et-st) << " ms" << endl;
    cout<<"Для продолжения работы программы нажмите любую клавишу...";

    cin.get();
// return rc;
}
