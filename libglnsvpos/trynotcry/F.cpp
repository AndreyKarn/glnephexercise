#include <math.h>
#include <iostream>
#include <stdio.h>


using namespace std;
//константы
const double J02=1082625.75e-9;
const double GM=398600441.8e6;
const double ae=6378136;

struct coord
{
    double xa, ya, za,vxa,vya,vza;
};

coord  F(struct coord cor, double xa ,double ya,double za,double vxa,double vya,double vza)
{


double r, GMrat, xarat, yarat, zarat, ro, Vxa, Vya, Vza ;
r=sqrt(pow(xa,2)+pow(ya,2)+pow(za,2));
GMrat=GM/(pow(r,2));
xarat=xa/r;
yarat=ya/r;
zarat=za/r;
ro=ae/r;

Vxa=-GMrat*xarat-(3/2)*J02*GMrat*xarat*(pow(ro,2))*(1-5*pow(zarat,2));//  %+Jxas+Jxam;
Vya=-GMrat*yarat-(3/2)*J02*GMrat*yarat*(pow(ro,2))*(1-5*pow(zarat,2));//%+Jyas+Jyam;
Vza=-GMrat*zarat-(3/2)*J02*GMrat*zarat*(pow(ro,2))*(3-5*pow(zarat,2));//%+Jzas+Jzam;
cor.xa=vxa;
cor.ya=vya;
cor.za=vza;
cor.vxa=Vxa;
cor.vya=Vya;
cor.vza=Vza;
return cor;
}

