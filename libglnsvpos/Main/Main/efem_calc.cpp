#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include"efem_calc.h"
#include"Optimization.h"


/*
������ ������� ������� ������������ ������ �������� � ��������� ��������� ��
� ������������ ��������������� ��, ���������� �� ��� ���������� ��������� �� ������� �� ������ �����-�����
�� �������� ������� ������� � ������������ ���������� ��������� �� �������
�� ��������� � ������������ �������������� � MatLab

������� ������:
double x = 10584969.2383; // ���������� � � ������� ��_90[m]
double y = 2721713.86719; // ���������� Y � ������� ��_90[m]
double z = 23027096.6797; // ���������� Z � ������� ��_90[m]
double V_x = -788.876533508; // �������� �� � � ������� ��_90[m / s]
double V_y = 3058.68911743; // �������� �� Y � ������� ��_90[m / s]
double V_z = 1.50871276855; // �������� �� Z � ������� ��_90[m / s]
double a_x = 0; // ��������� �� � � ������� ��_90[m / s ^ 2]
double a_y = 0; // ��������� �� Y � ������� ��_90[m / s ^ 2]
double a_z = -372.529029846e-9; // ��������� �� Z � ������� ��_90[m / s ^ 2]
int Year = 2020; // ��� �� UTC
int Month = 2; // ����� �� UTC
int Day = 10; // ���� �� UTC
int Hour = 13; // ��� �� UTC
int Min = 45; // ������ �� UTC
int Sec = 0; // ������� �� UTC
double t_start = 12; // ����� ������ �������� ��������� � ���
double t_end = 24; // ����� ��������� �������� ��������� � ���
int Day_start = 10; // ���� ������ ��������(������ ��������� � ���� ������� ��������) ��������� � ���
double delta_t = 0.1; // ��� �������

�������� ������:
���, �� ������� ��������������� ���������� ��������� �� �� �������� 
������ ������� � ������������ �� � MatLab � � �++ � ����� ���������� ���������

*/

using namespace std;


typedef struct Place
{
    double x;
    double y;
    double z;
    double V_x;
    double V_y;
    double V_z;

}Place;


void efem_calc(double x, double y, double z, double V_x, double V_y, double V_z,
    double a_x, double a_y, double a_z, int Year, int Month,
    int Day, int Hour, int Min, int Sec,
    int t_start, int t_end, double delta_t, int Day_start)

{
    if (Day_start != Day) {
        cout << "ERROR! The forecast day does not coincide with the day of arrival of the ephemeris\n\n";
        exit(0);
    }
    double omega_z = 7.2921151467e-5; // ������� ������� �������� �������� ����� ������������ �.��[��� / �]
    double A_e = 6378136; // ������� ������� ����������� ����������[�]
    double GM = 398600441.8e6; // ��������.�����.����.���� ����� � ������ ������.[� ^ 3 / c ^ 2]
    double J_2 = 1082625.75e-9; // ��������� ������������� ����.������ �������
    double pi = 3.141592653589793;
    // ������� � ��� ������ ��������
    t_start = (t_start + 3) * 60 * 60;// ����� ������ �������� ��������� � ���
    t_end = (t_end + 3) * 60 * 60; // ����� ��������� �������� ��������� � ���
    // �������� � ��� ������� ����������� ������
    if (Hour + 3 < 24) { // ���� �������� �� ���
        Hour = Hour + 3;
    }
    else {
        Hour = (Hour + 3) - 24;
        Day = Day + 1;
    }
    int  DEL = Year - 1995;
    int N_ch = ceil((DEL) / 4.0); // ����� ������������ �� ����.
    int N_t = 0;
    if (DEL % 4 == 1) {
        int vis[13] = { 31, 31 + 29, 31 + 29 + 31, 31 + 29 + 31, 31 + 29 + 31 + 30, 31 + 29 + 31 + 30 + 31,
            31 + 29 + 31 + 30 + 31 + 30, 31 + 29 + 31 + 30 + 31 + 30 + 31, 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31,
            31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30, 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31,
            31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30, 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + 31 };
        N_t = N_t + vis[Month - 2] + Day; // ����� ����� �� ������ ������������ �� ����.
    }
    else {
        N_t = N_t + 365 * ((DEL % 4) - 1) + 1;
        int n_vis[13] = { 31, 31 + 29, 31 + 28 + 31, 31 + 29 + 31, 31 + 29 + 31 + 30, 31 + 29 + 31 + 30 + 31,
            31 + 29 + 31 + 30 + 31 + 30, 31 + 29 + 31 + 30 + 31 + 30 + 31, 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31,
            31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30, 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31,
            31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30, 31 + 29 + 31 + 30 + 31 + 30 + 31 + 31 + 30 + 31 + 30 + 31 };
        N_t = N_t + n_vis[Month - 2] + Day; // ����� ����� �� ������ ������������ �� ����.
    }
    int t_e = Hour * 60 * 60 + Min * 60 + Sec; // ����� ������� �� ������ ����� �� ����.
    // ������������� ���������, ��������� � ��������� �� �� - 90 � �����.�������������� ��
    // �������� ������� �������� ����� �� ��������
    double JD0 = 1461 * ((double)N_ch - 1) + (double)N_t + 2450082.5; // ������� ��������� ���� �� 0 ����� ����� ���
    double ERA = 2 * pi * (0.7790572732640 + 1.00273781191135448 * (JD0 - 2451545)); //  ���� �������� �����[���]
    double delta_T = (JD0 - 2451545) / 36525; // ����� �� ����� 2000 � 1 ������ 12 �(UTC(SU))
    // �� ������� ����� � ��������� ��������� �� 36525 ����������� �����
    double GMST = ERA + 0.0000000703270726 + 0.0223603658710194 * delta_T \
        + 0.0000067465784654 * pow(delta_T, 2) - 0.0000000000021332 * pow(delta_T, 3) \
        - 0.0000000001452308 * pow(delta_T, 4) - 0.000000000001784 * pow(delta_T, 5); //  ������� �� - � ����� �� ��������[���]
    double S = GMST + omega_z * ((double)t_e - 10800);
    x = x * cos(S) - y * sin(S);
    y = x * sin(S) + y * cos(S);
    V_x = V_x * cos(S) - V_y * sin(S) - omega_z * y;
    V_y = V_x * sin(S) + V_y * cos(S) + omega_z * x;
    double JJ_x = a_x * cos(S) - a_y * sin(S);
    double JJ_y = a_x * sin(S) + a_y * cos(S);
    double JJ_z = a_z;



    int length = ((double)t_end - (double)t_start) / delta_t;
    Place* result = new Place[length];

    result[0].x = x;
    result[0].y = y;
    result[0].z = z;
    result[0].V_x = V_x;
    result[0].V_y = V_y;
    result[0].V_z = V_z;

    result = Optim(result, length, J_2,
        JJ_x, JJ_y, JJ_z,
        A_e, GM, delta_t, t_start, t_e, t_end); // ��������� �� ������� �� �������� 
    //������� ������� � ������������ ��������������� ��

    // ������ ������ �� .txt � ������� �������� ������ ��������� �� �� 
    //�������� ������� ������� � ������������ ��������������� �� ����������� � MatLab
    // � ��������� � ������������ � �++
    ifstream File;
    File.open("XYZ.txt");
    if (!File.is_open()) {
        cout << "Error opening file\n\n";
        exit(0);
    }
    else {
        cout << "File is opening\n\n";
    }

    double delta_X0, delta_Y0, delta_Z0, delta_X, delta_Y, delta_Z;
    double X_file, Y_file, Z_file;
    delta_X = 0;
    delta_Y = 0;
    delta_Z = 0;

    for (int i = 0; i < length; i++) {
        File >> X_file;
        delta_X0 = abs(X_file - result[i].x);
        File >> Y_file;
        delta_Y0 = abs(Y_file - result[i].y);
        File >> Z_file;
        delta_Z0 = abs(Z_file - result[i].z);
        if (delta_X0 > delta_X) {
            delta_X = delta_X0;
        }
        if (delta_Y0 > delta_Y) {
            delta_Y = delta_Y0;
        }
        if (delta_Z0 > delta_Z) {
            delta_Z = delta_Z0;
        }
    }
    std::cout.setf(std::ios::fixed);
    std::cout.precision(14); //14 - ����� �������� ����� �����
    std::cout << "Maximum calculation error in C++ compared to MatLab:\n\n";
    std::cout << "delta_X, m = " << delta_X << "\n\n";
    std::cout << "delta_Y, m = " << delta_Y << "\n\n";
    std::cout << "delta_Z, m = " << delta_Z << "\n\n";
    float end_time = clock() / 1000.0;
    std::cout.setf(std::ios::fixed);
    std::cout.precision(3); //3 - ����� �������� ����� �����
    std::cout << "Time of the programm:\n\n";
    std::cout << "Time, s = " << end_time << "\n\n";
}

