clear all; clc; close all;
format long g
%time
%
%���������
pi=3.14159265359;
% J02=1082625.75*10^-9;  %��������� ������������� ����������� ������ �������
 ae=6378136;         %������� (��������������) ������� ����������� ����������
% GM=398600441.8*10^6;    %� ��������������� ��������� ��������������� ���� ����� 
% J20=1082625.75*10^-9;    %��������� ������������� ����������� ������ �������, ��������������� �������� ������ �����
% we=0.7292115*10^-4;        %earth's rotation rate
we=0.7292115*10^-4; %earth's rotation rate
%������ �� RTKNAVI

x0=10192674.32;
y0=-12367565.43;
z0=19866879.39;
vx=2599.78676;
vy=-789.66141;
vz=-1827.75784;
ax=0.0000019;
ay=0.000009;
az=-0.0000028;
Tau=-38310; %ns
Gamma=0.0018; %ns
%%
%time format
 %2020.02.10 13.45.18           !!!!���������!!!!!
 N4=7;
 Nt=41;
 te=(13+3)*60*60+45*60;
 
% ������ � ����� ������� ����������
time_start=(12+3)*60*60; %(+3 UTC)
time_final=(24+3)*60*60;

C1=0;%use after 2119 year
C2=0;%use after 2239 year
JD0=1461*(N4-1)+Nt+2450082.5-(Nt-3)/25+C1+C2; %������� ��������� ���� �� 0 ����� ����� ���
T=(JD0+(te -10800)/86400-2451545.0)/36525;

%������� ������� ������ GMST � ������ (���������� � ���)
JDN=JD0+0.5;

Tdel=(JD0-2451545.0)/36525;
ERA=2*pi*(0.7790572732640+1.00273781191135448*(JD0-2451545.0));%���� �������� �����, ���
GMST=ERA+0.0000000703270726+0.0223603658710194*Tdel+0.0000067465784654*Tdel^2-0.0000000000021332*Tdel^3-0.0000000001452308*Tdel^4-0.0000000000001784*Tdel^5;   %�������� �������� ����� �� �������� (���) (GST ���)
S=GMST+we*(te-10800);  %10800 �� ���

%%
%Coordinates transformation to an inertial reference frame:

%Position
xa=x0*cos(S)-y0*sin(S);
ya=x0*sin(S)+y0*cos(S);
za=z0;
%Velocity
vxa=vx*cos(S)-vy*sin(S)-we*ya;
vya=vx*sin(S)+vy*cos(S)+we*xa;
vza=vz;

Jsm_x=ax*cos(S)-ay*sin(S);
Jsm_y=ax*sin(S)+ay*cos(S);
Jsm_z=az;

%%
%��������

coordinat=math(xa,ya,za,vxa,vya,vza,Jsm_x,Jsm_y,Jsm_z,time_start,time_final, te,T);

%������ �����
[EAR_x,EAR_y,EAR_z] = sphere(20);
EAR_x=ae.*EAR_x;
EAR_y=ae.*EAR_y;
EAR_z=ae.*EAR_z;
figure (1)


surf(EAR_x,EAR_y,EAR_z)
hold on
grid on
plot3(coordinat(:,1),coordinat(:,2),coordinat(:,3))