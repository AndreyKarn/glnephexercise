#include <iostream>
#include <vector>
#include <array>
#include <math.h>
const double pi=3.14159265359;
const double we=7.2921151467e-5; //earth's rotation rate


struct tm
{
    double S,time_start,time_final,T,te,GMST,ERA,JD0;
};

tm times(int N4, int Nt, int hour, int minut, int sec, int h_st, int h_fin, //вход
            struct tm time  )//выход
{
int Tdel;
te=(hour+3)*60*60+60*minut+sec;//+3UTC
time_start=(h_st+3)*60*60; //(+3 UTC)
time_final=(h_fin+3)*60*60;
JD0=1461*(N4-1)+Nt+2450082.5;//текущая юлианская дата на 0 часов шкалы МДВ
T=(JD0+(te -10800)/86400-2451545.0)/36525;
//расчеты времени всякие GMST и прочее (Приложение Л ИКД)


Tdel=(JD0-2451545.0)/36525;

ERA=2*pi*(0.7790572732640 + 1.00273781191135448*(JD0 - 2451545.0));//%угол поворота Земли, рад
GMST=ERA+0.0000000703270726+0.0223603658710194*Tdel+0.0000067465784654*Tdel*Tdel-0.0000000000021332*Tdel*Tdel*Tdel-0.0000000001452308*Tdel*Tdel*Tdel*Tdel-0.0000000000001784*Tdel*Tdel*Tdel*Tdel*Tdel;  // %истинное звездное время по Гринвичу (рад) (GST ИКД)
S=GMST+we*(te-10800);  //%10800 из ИКД
}
