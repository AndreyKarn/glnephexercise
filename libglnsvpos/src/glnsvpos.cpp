#include <include\libglnsvpos\glnsvpos.h>
#include <include\libglnsvpos\rungekutta.h>

#include <iostream>
#include <cmath>
#include <ostream>

using namespace std;

void koordinate_GLONASS(double **koord_n)
{
    double *koord_0 = new double [6];

    // начальные условия
    koord_0[0] = 2656202.15;
    koord_0[1] = 19596105.96;
    koord_0[2] = 16152160.15;
    koord_0[3] = -453.223228;
    koord_0[4] = 2162.4279;
    koord_0[5] = -2549.744606;

    double Ax = -0.0000019;
    double Ay = -0.0000019;
    double Az = -0.0000019;
    double w_e = 7.2921151467e-5; //средняя угловая скорость вращения Земли

    //Расчет времени формата ГЛОНАСС
    //Дата: 10.02.20 Время: 13:45:18
    double h = 0.1;
    double N4 = (2020-1996)/4.0 + 1; //Номер четырехлетнего периода
    double Nt = 31 + 25 + 1; //Текущие сутки от начала года
    double tb = 18 + 45*60 + 13*60*60 + 10800; //Момент времени по шкале МДВ
    double t_start = 12; //Время начала прогноза
    double t_end = 24; //Время окончания
    double T_start = (t_start + 3) * 60 * 60; //по шкале МДВ
    double T_end = (t_end + 3) * 60 * 60; //по шкале МДВ
    double JD0 = 1461 * (N4 - 1) + Nt + 2450082.5 - (Nt - 3)/25.0;
    double del_T = (JD0 - 2451545)/36525.0;
    double ERA = 2*M_PI * (0.7790572732640 + 1.00273781191135448 * (JD0 - 2451545));
    double GMST = ERA + 0.0000000703270726 + 0.0223603658710194 * del_T + 0.0000067465784654 * pow(del_T,2) - 0.0000000000021332 * pow(del_T,3) - 0.0000000001452308 * pow(del_T,4) - 0.0000000000001784 * pow(del_T,5);

    int num = (int) abs((T_end - T_start)/h)+1;
    double *ti = new double [num];
    for (int i = 0; i < num; i++)
    {
        ti[i] = T_start + i*h;
    }
    int num_eph = (int) (tb - T_start)/h;
    double S = GMST + w_e * (tb - 10800);
    double cos_S = cos(S);
    double sin_S = sin(S);
    // пересчет координат из ПЗ-90 в ECI
    koord_n[num_eph][0] = koord_0[0]*cos_S - koord_0[1]*sin_S;
    koord_n[num_eph][1] = koord_0[0]*sin_S + koord_0[1]*cos_S;
    koord_n[num_eph][2] = koord_0[2];
    koord_n[num_eph][3] = koord_0[3]*cos_S - koord_0[4]*sin_S - w_e*koord_n[num_eph][1];
    koord_n[num_eph][4] = koord_0[3]*sin_S + koord_0[4]*cos_S + w_e*koord_n[num_eph][0];
    koord_n[num_eph][5] = koord_0[5];
    double Ax0 = Ax * cos_S - Ay * sin_S;
    double Ay0 = Ax * sin_S + Ay * cos_S;
    double Az0 = Az;

    for (int i = num_eph; i > 0; i--)
    {
        RungeKutta(ti[i], ti[i-1], koord_n[i], koord_n[i-1]);
    }
    for (int i = num_eph; i < num-1 ; i++)
    {
        RungeKutta(ti[i], ti[i+1], koord_n[i], koord_n[i+1]);
    }
    for (int i = 0; i < num; i++)
    {
        double tau = ti[i] - tb;
        double del_X = Ax0 * tau * tau/2.0;
        double del_Y = Ay0 * tau * tau/2.0;
        double del_Z = Az0 * tau * tau/2.0;
        double del_Vx = Ax0 * tau;
        double del_Vy = Ay0 * tau;
        double del_Vz = Az0 * tau;
        //Добавляем поправки к результатам интегрирования
        koord_n[i][0] += del_X;
        koord_n[i][1] += del_Y;
        koord_n[i][2] += del_Z;
        koord_n[i][3] += del_Vx;
        koord_n[i][4] += del_Vy;
        koord_n[i][5] += del_Vz;
    }
    delete[] koord_0;
    koord_0 = nullptr;
    delete[] ti;
    ti = nullptr;
}
