#include <include/libglnsvpos/rungekutta.h>

using namespace std;

Y_s* diffs(double tn , struct Y_s Y)
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

    struct Y_s* dY;
    dY = new struct Y_s;
    // Дифуры
    dY->F1 = Y.F4;
    dY->F2 = Y.F5;
    dY->F3 = Y.F6;

    dY->F4 = - GM0 * crdX0 - 3/2 * J02 * GM0 * crdX0 * Rho * Rho * (1 - 5 * crdZ0 * crdZ0);
    dY->F5 = - GM0 * crdY0 - 3/2 * J02 * GM0 * crdY0 * Rho * Rho * (1 - 5 * crdZ0 * crdZ0);
    dY->F6 = - GM0 * crdZ0 - 3/2 * J02 * GM0 * crdZ0 * Rho * Rho * (3 - 5 * crdZ0 * crdZ0);

    return dY;
}

int RK(uint32_t N, double h, struct Y_s* Y) {

    if (N == 0 || h == 0) return 0;

    struct Y_s *k1, *k2, *k3, *k4, *knextstep;
    struct Y_s Y2, Y3, Y4;

    for (uint32_t k = 1; k < N; k++)
    {
        k1 = new struct Y_s;
        k2 = new struct Y_s;
        k3 = new struct Y_s;
        k4 = new struct Y_s;
        knextstep = new struct Y_s;

        k1 = diffs(0, Y[k-1]);

        Y2.F1 = Y->F1 + h * k1->F1 / 2;
        Y2.F2 = Y->F2 + h * k1->F2 / 2;
        Y2.F3 = Y->F3 + h * k1->F3 / 2;
        Y2.F4 = Y->F4 + h * k1->F4 / 2;
        Y2.F5 = Y->F5 + h * k1->F5 / 2;
        Y2.F6 = Y->F6 + h * k1->F6 / 2;

        k2 = diffs (0 + h / 2, Y2);

        Y3.F1 = Y->F1 + h * k2->F1 / 2;
        Y3.F2 = Y->F2 + h * k2->F2 / 2;
        Y3.F3 = Y->F3 + h * k2->F3 / 2;
        Y3.F4 = Y->F4 + h * k2->F4 / 2;
        Y3.F5 = Y->F5 + h * k2->F5 / 2;
        Y3.F6 = Y->F6 + h * k2->F6 / 2;

        k3 = diffs (0 + h / 2 , Y3);

        Y4.F1 = Y->F1 + h * k2->F1;
        Y4.F2 = Y->F2 + h * k2->F2;
        Y4.F3 = Y->F3 + h * k2->F3;
        Y4.F4 = Y->F4 + h * k2->F4;
        Y4.F5 = Y->F5 + h * k2->F5;
        Y4.F6 = Y->F6 + h * k2->F6;

        k4 = diffs (0 + h , Y4);

        knextstep->F1 = (double)h / 6 * ( k1->F1 + 2 * k2->F1 + 2 * k3->F1 + k4->F1 );
        knextstep->F2 = (double)h / 6 * ( k1->F2 + 2 * k2->F2 + 2 * k3->F2 + k4->F2 );
        knextstep->F3 = (double)h / 6 * ( k1->F3 + 2 * k2->F3 + 2 * k3->F3 + k4->F3 );
        knextstep->F4 = (double)h / 6 * ( k1->F4 + 2 * k2->F4 + 2 * k3->F4 + k4->F4 );
        knextstep->F5 = (double)h / 6 * ( k1->F5 + 2 * k2->F5 + 2 * k3->F5 + k4->F5 );
        knextstep->F6 = (double)h / 6 * ( k1->F6 + 2 * k2->F6 + 2 * k3->F6 + k4->F6 );

        Y[k].F1 = Y[k-1].F1 + knextstep->F1;
        Y[k].F2 = Y[k-1].F2 + knextstep->F2;
        Y[k].F3 = Y[k-1].F3 + knextstep->F3;
        Y[k].F4 = Y[k-1].F4 + knextstep->F4;
        Y[k].F5 = Y[k-1].F5 + knextstep->F5;
        Y[k].F6 = Y[k-1].F6 + knextstep->F6;

        cout << "Y[" << k <<"].F1 = " << Y[k].F1 << endl;

        delete k1;
        delete k2;
        delete k3;
        delete k4;
        delete knextstep;
    }
    return 0;
};

int mult(int a, int b){
    return a*b;
}
