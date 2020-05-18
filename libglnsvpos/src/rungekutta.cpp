#include <include/libglnsvpos/rungekutta.h>

using namespace std;

Y_s* diffs(double tn , struct Y_s Y)
{
    double J02 = 1082625.75e-9; // зональный гармонический коэффициент второй степени, характеризующий полярное сжатие Земли
    double GM = 398600441.8e6; // геоцентрическая константа гравитационного поля Земли с учетом атмосферы, [м3/c2]
    double a_e = 6378136; // большая полуось общеземного эллипсоида, [м]

    double crdX = Y.X;
    double crdY = Y.Y;
    double crdZ = Y.Z;

    double r = sqrt(crdX * crdX + crdY * crdY + crdZ * crdZ);

    double GM0 = GM / (r * r);
    double Rho = a_e / r;
    double crdX0 = crdX / r;
    double crdY0 = crdY / r;
    double crdZ0 = crdZ / r;

    struct Y_s* dY;
    dY = new struct Y_s;
    // Дифуры
    dY->X = Y.VX;
    dY->Y = Y.VY;
    dY->Z = Y.VZ;

    dY->VX = - GM0 * crdX0 - (double)1.5 * J02 * GM0 * crdX0 * Rho * Rho * (1 - (double)5 * crdZ0 * crdZ0);
    dY->VY = - GM0 * crdY0 - (double)1.5 * J02 * GM0 * crdY0 * Rho * Rho * (1 - (double)5 * crdZ0 * crdZ0);
    dY->VZ = - GM0 * crdZ0 - (double)1.5 * J02 * GM0 * crdZ0 * Rho * Rho * (3 - (double)5 * crdZ0 * crdZ0);

    return dY;
}

int RK(uint32_t N, double h, struct Y_s* Y) {

    if (N == 0 || h == 0) return 0;

    struct Y_s *k1, *k2, *k3, *k4, *knextstep;
    struct Y_s Y2, Y3, Y4;

    for (uint32_t k = 1; k <= N; k++)
    {
        k1 = new struct Y_s;
        k2 = new struct Y_s;
        k3 = new struct Y_s;
        k4 = new struct Y_s;
        knextstep = new struct Y_s;

        k1 = diffs(0, Y[k-1]);

        Y2.X = Y[k-1].X + h * k1->X / (double)2;
        Y2.Y = Y[k-1].Y + h * k1->Y / (double)2;
        Y2.Z = Y[k-1].Z + h * k1->Z / (double)2;
        Y2.VX = Y[k-1].VX + h * k1->VX / (double)2;
        Y2.VY = Y[k-1].VY + h * k1->VY / (double)2;
        Y2.VZ = Y[k-1].VZ + h * k1->VZ / (double)2;

        k2 = diffs (0 + h / 2, Y2);

        Y3.X = Y[k-1].X + h * k2->X / (double)2;
        Y3.Y = Y[k-1].Y + h * k2->Y / (double)2;
        Y3.Z = Y[k-1].Z + h * k2->Z / (double)2;
        Y3.VX = Y[k-1].VX + h * k2->VX / (double)2;
        Y3.VY = Y[k-1].VY + h * k2->VY / (double)2;
        Y3.VZ = Y[k-1].VZ + h * k2->VZ / (double)2;

        k3 = diffs (0 + h / 2 , Y3);

        Y4.X = Y[k-1].X + h * k3->X;
        Y4.Y = Y[k-1].Y + h * k3->Y;
        Y4.Z = Y[k-1].Z + h * k3->Z;
        Y4.VX = Y[k-1].VX + h * k3->VX;
        Y4.VY = Y[k-1].VY + h * k3->VY;
        Y4.VZ = Y[k-1].VZ + h * k3->VZ;

        k4 = diffs (0 + h , Y4);

        knextstep->X = h / (double)6 * ( k1->X + 2 * k2->X + 2 * k3->X + k4->X );
        knextstep->Y = h / (double)6 * ( k1->Y + 2 * k2->Y + 2 * k3->Y + k4->Y );
        knextstep->Z = h / (double)6 * ( k1->Z + 2 * k2->Z + 2 * k3->Z + k4->Z );
        knextstep->VX = h / (double)6 * ( k1->VX + 2 * k2->VX + 2 * k3->VX + k4->VX );
        knextstep->VY = h / (double)6 * ( k1->VY + 2 * k2->VY + 2 * k3->VY + k4->VY );
        knextstep->VZ = h / (double)6 * ( k1->VZ + 2 * k2->VZ + 2 * k3->VZ + k4->VZ );

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

int mult(int a, int b){
    return a*b;
}
