clear all; clc;
format long;

%�������� ������

%�������� ���������:
%����������
X = 2656202.15;
Y = 19596105.96;
Z = 16152160.15;
%���������� ������� ��������
Vx = -453.223228;
Vy = 2162.4279;
Vz = -2549.744606;
%�����-��������� ���������
Ax = -0.0000019;
Ay = -0.0000019;
Az = -0.0000019;

w_e = 7.2921151467e-5; %������� ������� �������� �������� �����


%������ ������� ������� �������
%����: 10.02.20 �����: 13:45:18
N4 = (2020-1996)/4 + 1; %����� �������������� �������
Nt = 31 + 25 + 1; %������� ����� �� ������ ����
tb = 18 + 45*60 + 13*60*60 + 10800; %������ ������� �� ����� ���
t_start = 12; %����� ������ ��������
t_end = 24; %����� ���������
T_start = (t_start + 3) * 60 * 60; %�� ����� ���
T_end = (t_end + 3) * 60 * 60; %�� ����� ���


%������ �������� ��������� ������� �� ��������
JD0 = 1461 * (N4 - 1) + Nt + 2450082.5 - (Nt - 3)/25;
del_T = (JD0 - 2451545)/36525;
ERA = 2*pi * (0.7790572732640 + 1.00273781191135448 * (JD0 - 2451545));
GMST = ERA + 0.0000000703270726 + 0.0223603658710194 * del_T + 0.0000067465784654 * del_T^2 - 0.0000000000021332 * del_T^3 - 0.0000000001452308 * del_T^4 - 0.0000000000001784 * del_T^5;


%�������� � ������������ ��������������� ������� ���������
S = GMST + w_e * (tb - 10800);

X0 = X * cos(S) - Y * sin(S);
Y0 = X * sin(S) + Y * cos(S);
Z0 = Z;

Vx0 = Vx * cos(S) - Vy * sin(S) - w_e * Y0;
Vy0 = Vx * sin(S) + Vy * cos(S) + w_e * X0;
Vz0 = Vz;


%�������������� ��������� �������
h = 1; %������ ��� 
ti = T_start:h:T_end; 
RK = [X0 Y0 Z0 Vx0 Vy0 Vz0]; %������ ��������� ��������� �������

[t, f] = ode45('integ', tb:-h:ti(1), RK); %����� �����-����� 4 �������
tt = t(end:-1:2,:); 
ff = f(end:-1:2,:);
[t, f] = ode45('integ', tb:h:ti(end), RK);
tt = [tt;t]; %������� ������ ������������ �����
ff = [ff; f]; %������� � ���������



%���������� ��������-������ ���������
Ax0 = Ax * cos(S) - Ay * sin(S);
Ay0 = Ax * sin(S) + Ay * cos(S);
Az0 = Az;

tau = tt - tb;
del_X = Ax0 * (tau.^2)/2;
del_Y = Ay0 * (tau.^2)/2;
del_Z = Az0 * (tau.^2)/2;
del_Vx = Ax0 * tau;
del_Vy = Ay0 * tau;
del_Vz = Az0 * tau;
del = [del_X del_Y del_Z del_Vx del_Vy del_Vz];


%��������� �������� � ����������� ��������������
ff = ff + del; 


%�������� ��������� � ��������� � ������ ������� ��������� ��-90
SS = GMST + w_e * (tt - 10800);
coor(:,1) = ff(:,1) .* cos(SS) + ff(:,2) .* sin(SS);
coor(:,2) = -ff(:,1) .* sin(SS) + ff(:,2) .* cos(SS);
coor(:,3) = ff(:,3);


%�������� ��������� ������� � � ������� WGS-84
%������� - 55 �������� 45 ����� 24 �������
%������ - 37 �������� 42 ������ 11 ������
L = 0.97313; %������� � ��������
H = 0.65804; %������ � ��������
E = 500; %������
PRM = [L H E];


%���������� SkyPlot
for i = 1:length(coor(:,1))

 [x(i) y(i) z(i)] = ecef2enu(coor(i,1),coor(i,2),coor(i,3),L,H,E,wgs84Ellipsoid,'radians');
 
 if z(i) > 0
 r(i) = sqrt(x(i)^2 + y(i)^2 + z(i)^2);
 teta(i) = acos(z(i)/r(i));

 if x(i) > 0
     phi(i) = -atan(y(i)/x(i))+pi/2;
 elseif (x(i)<0)&&(y(i)>0)
 phi(i) = -atan(y(i)/x(i))+3*pi/2;
 elseif (x(i)<0)&&(y(i)<0)
 phi(i) = -atan(y(i)/x(i))-pi/2;
 end
 
 else teta(i) = NaN;
 r(i) = NaN;
 phi(i) = NaN;

 end
end


%�������
Rz = 6371e3;
[Xz,Yz,Zz] = sphere(30);
figure(1)
surf(Rz*Xz,Rz*Yz,Rz*Zz)
grid on
hold on
plot3(ff(:,1),ff(:,2),ff(:,3),'g')
plot3(coor(:,1),coor(:,2),coor(:,3),'b')
xlabel('�,�')
ylabel('Y,�')
zlabel('Z,�')


figure(2)
ax = polaraxes;
polarplot(ax,phi,teta*180/pi)
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
title('SkyPlot �� �12')








