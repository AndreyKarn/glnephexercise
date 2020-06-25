#include <iostream>
#include <tuple>
#include <fstream>
#include <ctime>
#include <cmath>
#include <rk.h>
#include <integ.h>


using namespace std;

int main()
{
    time_t start, end;
    //Исходные данные
    //Заданные эфемериды:
    //Координаты
    double X = 2656202.15;
    double Y = 19596105.96;
    double Z = 16152160.15;
    //Компоненты вектора скорости
    double Vx = -453.223228;
    double Vy = 2162.4279;
    double Vz = -2549.744606;
    //Лунно-солнечные ускорения
    double Ax = -0.0000019;
    double Ay = -0.0000019;
    double Az = -0.0000019;

    double w_e = 7.2921151467e-5; //средняя угловая скорость вращения Земли


    //Расчет времени формата ГЛОНАСС
    //Дата: 10.02.20 Время: 13:45:18
    int N4 = (2020-1996)/4.0 + 1; //Номер четырехлетнего периода
    int Nt = 31 + 25 + 1; //Текущие сутки от начала года
    int tb = 18 + 45*60 + 13*60*60 + 10800; //Момент времени по шкале МДВ
    int T_start = (12 + 3) * 60 * 60; //по шкале МДВ
    int T_end = (24 + 3) * 60 * 60; //по шкале МДВ


    //Расчет среднего звездного времени по Гринвичу
    double JD0 = 1461 * (N4 - 1) + Nt + 2450082.5 - (Nt - 3.0)/25.0;
    double del_T = (JD0 - 2451545)/36525.0;
    double ERA = 2*M_PI * (0.7790572732640 + 1.00273781191135448 * (JD0 - 2451545));
    double GMST = ERA + 0.0000000703270726 + 0.0223603658710194 * del_T + 0.0000067465784654 * pow(del_T,2) - 0.0000000000021332 * pow(del_T,3) - 0.0000000001452308 * pow(del_T,4) - 0.0000000000001784 * pow(del_T,5);


    //Пересчет в инерциальную геоцентрическую систему координат
    double S = GMST + w_e * (tb - 10800);
    double sin_S = sin(S);
    double cos_S = cos(S);

    double X0 = X * cos_S - Y * sin_S;
    double Y0 = X * sin_S + Y * cos_S;
    double Z0 = Z;

    double Vx0 = Vx * cos_S - Vy * sin_S - w_e * Y0;
    double Vy0 = Vx * sin_S + Vy * cos_S + w_e * X0;
    double Vz0 = Vz;

    double Ax0 = Ax * cos_S - Ay * sin_S;
    double Ay0 = Ax * sin_S + Ay * cos_S;
    double Az0 = Az;

    //Интегрирование численным методом
    double h = 0.1; //задаем шаг
    int num = (int) (T_end - T_start)/h;
    int num_eph = (int) (tb - T_start)/h+1;
    double *ti = new double [num+1];
    double **ff = new double *[num+1];
    for (int i = 0; i < num+1; i++){
        ti[i] = T_start + i*h;
        ff[i] = new double [6];
    }
     //Вектор начальных состояний системы
    double *RK = new double [6] {X0, Y0, Z0, Vx0, Vy0, Vz0};
    for (int i = 0; i < 6; i++){
        ff[num_eph][i] = RK[i];
    }
    time(&start);
    for (int i = num_eph-1; i >= 0; i--){
        auto yn = RK4(ti[i], ti[i+1], ff[i+1]);
        ff[i][0] = get<0>(yn);
        ff[i][1] = get<1>(yn);
        ff[i][2] = get<2>(yn);
        ff[i][3] = get<3>(yn);
        ff[i][4] = get<4>(yn);
        ff[i][5] = get<5>(yn);
    }
    for (int i = num_eph+1; i < num+1; i++){
        auto yn = RK4(ti[i], ti[i-1], ff[i-1]);
        ff[i][0] = get<0>(yn);
        ff[i][1] = get<1>(yn);
        ff[i][2] = get<2>(yn);
        ff[i][3] = get<3>(yn);
        ff[i][4] = get<4>(yn);
        ff[i][5] = get<5>(yn);
    }
    for (int i = num_eph+1; i < num+1; i++){
        double tau = ti[i] - tb;
        ff[i][0] += Ax0 * tau*tau*0.5;
        ff[i][1] += Ay0 * tau*tau*0.5;
        ff[i][2] += Az0 * tau*tau*0.5;
        ff[i][3] += Ax0 * tau;
        ff[i][4] += Ay0 * tau;
        ff[i][5] += Az0 * tau;
    }
    time(&end);
    ofstream out;
    out.open("D:\\cpp_results.txt");
    ifstream in("D:\\matlab_results.txt");
    if (!in)
    {
        cout << "File not open, check file!" << endl;
    } else {
        cout << "File opened!" << endl;
    }
    double * ff_matlab = new double [3];
    double max_del = 0;
    int i_max = 0;
    for (int i = 0; i < num+1; i++)
    {
        in >> ff_matlab[0] >> ff_matlab[1] >> ff_matlab[2];
        string koord_str1 = to_string(ff[i][0]);
        string koord_str2 = to_string(ff[i][1]);
        string koord_str3 = to_string(ff[i][2]);
        out << koord_str1 << "\t" << koord_str2 << "\t" << koord_str3 << endl;
        for (int j = 0; j < 3; j++)
        {
            if (abs(ff[i][j] - ff_matlab[j]) > max_del)
            {
                max_del = abs(ff[i][j] - ff_matlab[j]);
                i_max = i;
            }
        }
    }
    time(&end);
    in.close();
    out.close();
    delete [] ti;
    ti = nullptr;
    delete [] RK;
    RK = nullptr;
    for (int i = 0; i < num+1; i++)
    {
        delete [] ff[i];
        ff[i] = nullptr;
    }
    delete [] ff;
    ff = nullptr;
    delete [] ff_matlab;
    ff_matlab = nullptr;
    double time_RK = difftime(end, start);
    string time_RK1 = to_string(time_RK*1000000/num);
    cout << "time of work of one cicle~: " << time_RK1 << " mcs" << endl;
    string max_del1 = to_string(max_del);
    cout << "max delta of coords: " << max_del1 << " m" << endl;
    string imax = to_string(i_max);
    cout << "number of max delta: " << imax << endl;
    return 0;
}
