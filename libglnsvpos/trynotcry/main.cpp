#include <iostream>
#include <vector>
#include <array>
#include <math.h>
#include <fstream>
//#include <golova.h>
const double pi=3.14159265359;
const double we=7.2921151467e-5; //earth's rotation rate

using namespace std;
//��������� �������
struct coord
{
    double xa, ya, za,vxa,vya,vza;
};




void RungKUTT(coord res[], double t, double dt);


void math_2(coord result[],int delt, int time_start,int time_final,int te,
            double xa, double ya,double za,double vxa, double vya, double vza, double Jsm_x, double Jsm_y,double Jsm_z);


int main()
{
double S,T,GMST,ERA,JD0, Tdel;
double N4,Nt,hour,minut,sec, h_st,h_fin ;
int time_start,time_final, te;
//��� ��� �����
coord *result;

N4=7;
 Nt=41;
 hour=13;
 minut=45;
 sec=0;
 h_st=12;  // %������ ����������, �����
 h_fin=24;  //%����� ����������, �����


te=(hour+3)*60*60+60*minut+sec;//+3UTC
time_start=(h_st+3)*60*60; //(+3 UTC)
time_final=(h_fin+3)*60*60;
JD0=1461*(N4-1)+Nt+2450082.5;//������� ��������� ���� �� 0 ����� ����� ���
T=(JD0+(te -10800)/86400-2451545.0)/36525;
//������� ������� ������ GMST � ������ (���������� � ���)
double TE=te;
Tdel=(JD0-2451545.0)/36525;
ERA=2*pi*(0.7790572732640 + 1.00273781191135448*(JD0 - 2451545.0));//%���� �������� �����, ���
GMST=ERA+0.0000000703270726+0.0223603658710194*Tdel+0.0000067465784654*Tdel*Tdel-0.0000000000021332*Tdel*Tdel*Tdel-0.0000000001452308*Tdel*Tdel*Tdel*Tdel-0.0000000000001784*Tdel*Tdel*Tdel*Tdel*Tdel;  // %�������� �������� ����� �� �������� (���) (GST ���)
S=GMST+we*(TE-10800);  //%10800 �� ���


// ��� ������ ������ ��������, ������������� ���

 double x0, y0,z0,vx,vy,vz,ax,ay,az,Tau,Gamma;
x0=10192674.32;
y0=-12367565.43;
z0=19866879.39;
vx=2599.78676;
vy=-789.66141;
vz=-1827.75784;
ax=0.0000019;
ay=0.000009;
az=-0.0000028;
Tau=-38310; //%ns
Gamma=0.0018; //%ns

// ��� ������������� ���������� � ��.������

double xa, ya, za, vxa,vya,vza,Jsm_x,Jsm_y,Jsm_z;
xa=x0*cos(S)-y0*sin(S);
ya=xa*sin(S)+y0*cos(S);
za=z0;
//%Velocity
vxa=vx*cos(S)-vy*sin(S)-we*ya;
vya=vx*sin(S)+vy*cos(S)+we*xa;
vza=vz;

Jsm_x=ax*cos(S)-ay*sin(S);
Jsm_y=ax*sin(S)+ay*cos(S);
Jsm_z=az;

int delt=43200;

result = new coord [delt];

math_2(result,delt,  time_start, time_final, te,
             xa,  ya, za, vxa, vya,  vza, Jsm_x, Jsm_y, Jsm_z);

ofstream output_x("res_x.txt");
for (int i = 0; i < delt; ++i)

{
output_x << result[i].xa <<'\n';
}


output_x.close();

ofstream output_y("res_y.txt");
for (int i = 0; i < delt; ++i)

{
output_y << result[i].ya <<'\n';
}


output_y.close();
ofstream output_z("res_z.txt");
for (int i = 0; i < delt; ++i)

{
output_z << result[i].za <<'\n';
}


output_z.close();


    std::cout <<"result[0].vxa   "  << result[0].vxa <<"\n";
    std::cout <<"result[1].vxa   "  << result[1].vxa <<"\n";
    std::cout <<"result[2].vxa   "  << result[2].vxa <<"\n";
    std::cout <<"result[3].vxa   "  << result[3].vxa <<"\n";

    delete []result;
}
