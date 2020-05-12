#include <include/libglnsvpos/rungekutta.h>
#include <stdio.h>

using namespace std;

Y_s diffs(double tn , struct Y_s Y)
{
    double J02 = 1082625.75e-9; // зональный гармонический коэффициент второй степени, характеризующий полярное сжатие Земли
    double GM = 398600441.8e6; // геоцентрическая константа гравитационного поля Земли с учетом атмосферы, [м3/c2]
    double a_e = 6378136; // большая полуось общеземного эллипсоида, [м]

    double crdX = Y.F1;
    double crdY = Y.F2;
    double crdZ = Y.F3;

    double r = sqrt(crdX * crdX + crdY * crdY + crdZ * crdZ);

    double GM0 = GM / (r * r);
    double Rho = a_e / r;
    double crdX0 = crdX / r;
    double crdY0 = crdY / r;
    double crdZ0 = crdZ / r;

    struct Y_s dY;
    // Дифуры
    dY.F1 = Y.F4;
    dY.F2 = Y.F5;
    dY.F3 = Y.F6;

    dY.F4 = - GM0 * crdX0 - 3/2 * J02 * GM0 * crdX0 * Rho * Rho * (1 - 5 * crdZ0 * crdZ0);
    dY.F5 = - GM0 * crdY0 - 3/2 * J02 * GM0 * crdY0 * Rho * Rho * (1 - 5 * crdZ0 * crdZ0);
    dY.F6 = - GM0 * crdZ0 - 3/2 * J02 * GM0 * crdZ0 * Rho * Rho * (3 - 5 * crdZ0 * crdZ0);

    return dY;
}

int RK(struct Ephemeris_s Eph) {

    // Коэф. фильтра
//    F0* k1 = new F0[N];
//    F0* k2 = new F0[N];
//    F0* k3 = new F0[N];
//    F0* k4 = new F0[N];

//    F0* Y = new F0[N];
//    F0* Y2 = new F0[N];
//    F0* Y3 = new F0[N];
//    F0* Y4 = new F0[N];


    // Начальные условия
//    Y->F1 = Eph.X;
//    Y->F2 = Eph.Y;
//    Y->F3 = Eph.Z;

//    Y->F4 = Eph.VX;
//    Y->F5 = Eph.VY;
//    Y->F6 = Eph.VZ;

//    cout << Y->F1 << std::endl; // value
//    cout << &(Y->F1) << std::endl; // adress

    //while (tn <= toe)
//    int k = 0;
//    while (k <= N)
//    {
//        k1 = diffs(tn, Y);

//        Y2.F1 = Y.F1 + h * k1.F1 / 2;
//        Y2.F2 = Y.F2 + h * k1.F2 / 2;
//        Y2.F3 = Y.F3 + h * k1.F3 / 2;
//        Y2.F4 = Y.F4 + h * k1.F4 / 2;
//        Y2.F5 = Y.F5 + h * k1.F5 / 2;
//        Y2.F6 = Y.F6 + h * k1.F6 / 2;

//        k2 = diffs (tn + h / 2 , Y2);

//        Y3.F1 = Y.F1 + h * k2.F1 / 2;
//        Y3.F2 = Y.F2 + h * k2.F2 / 2;
//        Y3.F3 = Y.F3 + h * k2.F3 / 2;
//        Y3.F4 = Y.F4 + h * k2.F4 / 2;
//        Y3.F5 = Y.F5 + h * k2.F5 / 2;
//        Y3.F6 = Y.F6 + h * k2.F6 / 2;

//        k3 = diffs (tn + h / 2 , Y3);

//        Y4.F1 = Y.F1 + h * k2.F1;
//        Y4.F2 = Y.F2 + h * k2.F2;
//        Y4.F3 = Y.F3 + h * k2.F3;
//        Y4.F4 = Y.F4 + h * k2.F4;
//        Y4.F5 = Y.F5 + h * k2.F5;
//        Y4.F6 = Y.F6 + h * k2.F6;

//        k4 = diffs (tn + h , Y4);

//        Y.F1 += h / 6 * ( k1.F1 + 2 * k2.F1 + 2 * k3.F1 + k4.F1 );
//        Y.F2 += h / 6 * ( k1.F2 + 2 * k2.F2 + 2 * k3.F2 + k4.F2 );
//        Y.F3 += h / 6 * ( k1.F3 + 2 * k2.F3 + 2 * k3.F3 + k4.F3 );
//        Y.F4 += h / 6 * ( k1.F4 + 2 * k2.F4 + 2 * k3.F4 + k4.F4 );
//        Y.F5 += h / 6 * ( k1.F5 + 2 * k2.F5 + 2 * k3.F5 + k4.F5 );
//        Y.F6 += h / 6 * ( k1.F6 + 2 * k2.F6 + 2 * k3.F6 + k4.F6 );
//        tn += h;
//    }
    return 0;
};

int mult(int a, int b){
    return a*b;
}
