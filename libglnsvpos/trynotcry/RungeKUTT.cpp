#include <math.h>
#include <iostream>
#include <iostream>
#include <vector>
#include <array>
#include <math.h>

struct coord
{
    double xa, ya, za,vxa,vya,vza;
};

coord  F(struct coord cor, double xa ,double ya,double za,double vxa,double vya,double vza);

//t- длина всего безобразия. dt- шаг. tst- начальное время для интегрирования
void RungKUTT(coord res[], double t, double dt)
{

coord K0;
coord K1;
coord K2;
coord K3;
coord ink1, ink2, ink3, ink;
coord cyk;

ink=res[0];
for (int i=1; i<t; i++)
{
    cyk=ink;

    K0=F(K0,cyk.xa,cyk.ya,cyk.za,cyk.vxa,cyk.vya,cyk.vza );

    ink1={cyk.xa+0.5*dt*K0.xa,cyk.ya+0.5*dt*K0.ya,cyk.za+0.5*dt*K0.za,cyk.vxa+0.5*dt*K0.vxa,cyk.vya+0.5*dt*K0.vya,cyk.vza+0.5*dt*K0.vza};       //(result(i,:)+0.5*dt.*K0));
    K1=F(K1, ink1.xa,ink1.ya,ink1.za,ink1.vxa,ink1.vya,ink1.vza );

    ink2={cyk.xa+0.5*dt*K1.xa,cyk.ya+0.5*dt*K1.ya,cyk.za+0.5*dt*K1.za,cyk.vxa+0.5*dt*K1.vxa,cyk.vya+0.5*dt*K1.vya,cyk.vza+0.5*dt*K1.vza};
    K2=F(K2, ink2.xa,ink2.ya,ink2.za,ink2.vxa,ink2.vya,ink2.vza);

    ink3={(cyk.xa+dt*K2.xa),(cyk.ya+dt*K2.ya),(cyk.za+dt*K2.za),(cyk.vxa+dt*K2.vxa),(cyk.vya+dt*K2.vya),(cyk.vza+dt*K2.vza)};
    K3=F(K3, ink3.xa,ink3.ya,ink3.za,ink3.vxa,ink3.vya,ink3.vza);

    ink={cyk.xa+(dt/6)*(K0.xa+2*K1.xa+2*K2.xa+K3.xa),cyk.ya+(dt/6)*(K0.ya+2*K1.ya+2*K2.ya+K3.ya),cyk.za+(dt/6)*(K0.za+2*K1.za+2*K2.za+K3.za),cyk.vxa+(dt/6)*(K0.vxa+2*K1.vxa+2*K2.vxa+K3.vxa),cyk.vya+(dt/6)*(K0.vya+2*K1.vya+2*K2.vya+K3.vya),cyk.vza+(dt/6)*(K0.vza+2*K1.vza+2*K2.vza+K3.vza)};

res[i]=ink;


}

}
