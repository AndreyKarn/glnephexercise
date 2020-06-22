#include <cmath>
#include <include\libglnsvpos\diffs.h>

void dif(double *koord, double *diffs)
{
    double J_2 = 1082625.75e-9; //зональный гармонический коэффициент второй степени, характеризующий полярное сжатие Земли
    double GM = 398600.4418e9; //геоцентрическая константа гравитационного поля Земли с учетом атмосферы
    double a_e = 6378136; //большая полуось общеземного эллипсоида
    double r = sqrt(koord[0]*koord[0] + koord[1]*koord[1] + koord[2]*koord[2]);
    diffs[0] = koord[3];
    diffs[1] = koord[4];
    diffs[2] = koord[5];
    double GM0 = GM/(r*r);
    double x0 = koord[0]/r;
    double y0 = koord[1]/r;
    double z0 = koord[2]/r;
    double R = a_e/r;
    diffs[3] = -GM0*x0 - 1.5*J_2*GM0*x0*R*R*(1-5*z0*z0);
    diffs[4] = -GM0*y0 - 1.5*J_2*GM0*y0*R*R*(1-5*z0*z0);
    diffs[5] = -GM0*z0 - 1.5*J_2*GM0*z0*R*R*(3-5*z0*z0);
}
