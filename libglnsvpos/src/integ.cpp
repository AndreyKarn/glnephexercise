#include<cmath>
#include<tuple>
#include<integ.h>


using namespace std;

//void integ(double *f, double *dif)
tuple<double,double,double,double,double,double> integ(double *f)
{
    double dX = f[3];
    double dY = f[4];
    double dZ = f[5];
    double J_2 = -1082625.75e-9;    //зональный гармонический коэффициент второй степени, характеризующий полярное сжатие Земли
    double GM = 398600.4418e9;    //геоцентрическая константа гравитационного поля Земли с учетом атмосферы
    double a_e = 6378136.0;         //большая полуось общеземного эллипсоида
    double r = sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]);
    double xx = f[0]/r;
    double yy = f[1]/r;
    double zz = f[2]/r;
    double Ro = a_e/r;
    double GM0 = GM/(r*r);
    double dVx = -GM0*xx - (3.0/2)*J_2*GM0*xx*Ro*Ro*(1-5*zz*zz);
    double dVy = -GM0*yy - (3.0/2)*J_2*GM0*yy*Ro*Ro*(1-5*zz*zz);
    double dVz = -GM0*zz - (3.0/2)*J_2*GM0*zz*Ro*Ro*(3-5*zz*zz);
    return make_tuple(dX,dY,dZ,dVx,dVy,dVz);
}
