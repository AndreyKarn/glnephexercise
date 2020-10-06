clear all, close all
tic
%внесем все константы
o_z = 7.2921151467e-5;
Rz=6371000;
J02 = -1082.63*10^-6;
ae = 6378136; 
GM = 398600441.8e6; 
N = 55.756687916667*pi/180;
E=37.703077194444*pi/180;
High = 500; 
%известные данные спутника номер 3
%координаты
Xkoor = 23036950.68; % м
Ykoor = -9091173.34; % м
Zkoor = 6041059.08; % м
%скорости
VXkoor= 755.86033; % м/с
VYkoor = -358.52718; % м/с
VZkoor = -3447.90649; % м/с
%ускорения
AXkoor = 0.0000056; % м/с^2
AYkoor = 0.0000000; % м/с^2
AZkoor = -0.0000028; % м/с^2
%время
%номер текущего четырёхлетия
N4=(2020-1996)/4+1;
%номер текущих суток
Nt=365*(2020-1996-4*(N4-1))+31+10+1;
%время  МДВ
time=13*60*60 + 45*60 + 18 + 10800;
%расчёт текущей юлианской даты
JD0 = 1461*(((2020-1996)/4)+1 - 1) + Nt + 2450008.5 - (Nt -3)/25; 
%среднее звёздное время по Гринвичу
T_d=(JD0-2451545)/36525;
ERA=2*pi*(0.7790572732640 + 1.00273781191135448*(JD0-2451545));
GMST=ERA+0.0000000703270726+0.0223603658710194*T_d+...
    +0.0000067465784654*T_d^2-0.0000000000021332*T_d^3+...
    - 0.0000000001452308*T_d^4-0.0000000000001784*T_d^5;
%пересчет в инерциальную систему координат 
S = GMST + o_z*(time - 3*60*60);
Xa=Xkoor*cos(S)-Ykoor*sin(S);
%координаты
Ya=Xkoor*sin(S)+Ykoor*cos(S);
Za=Zkoor;
%скорости
Vxa=VXkoor*cos(S)-VYkoor*sin(S)-o_z*Ya;
Vya=VXkoor*sin(S)+VYkoor*cos(S)+o_z*Xa;
Vza=VZkoor; 
%ускорения
Ax=AXkoor*cos(S)-AYkoor*sin(S);
Ay=AXkoor*sin(S)+AYkoor*cos(S);
Az=AZkoor;
%расчет 
r=sqrt(Xa^2 + Ya^2 + Za^2);
GMn = GM/r^2;
x1 = Xa/r;
y1 = Ya/r;
z1 = Za/r;
RO = ae/r;
%метод Рунге-Кутты
%задание времени и шага
TimeS = (12+3)*60*60;
TimeE = (24+3)*60*60;
del = 1;
pros = TimeS:del:TimeE;
F0 = [Xa Ya Za Vxa Vya Vza];
[t, F] = ode45('diffs', time:-del:pros(1), F0);
F1 = F(end:-1:2,:);
t1 = t(end:-1:2,:);
[t, F] = ode45('diffs', time:del:pros(end), F0);
F1 = [F1; F];
t1 = [t1;t];
% Учёт ускорений
tR1 = t1 - time;
AXR = AXkoor*(tR1.^2)/2;
AYR = AYkoor*(tR1.^2)/2;
AZR = AZkoor*(tR1.^2)/2; 
dVX = AXkoor*tR1;
dVY = AYkoor*tR1;
dVZ = AZkoor*tR1;
d_A = [AXR AYR AZR dVX dVY dVZ];
F1 = F1 + d_A;
% Пересчёт координат в ПЗ-90
S = GMST + o_z*(t1 - 3*60*60);
PZ90(:,1) = F1(:,1).*cos(S) + F1(:,2).*sin(S);
PZ90(:,2) = -F1(:,1).*sin(S) + F1(:,2).*cos(S);
PZ90(:,3) = F1(:,3);
cord_E = [N E High];
for i = 1:length(PZ90(:,1))
    [Xkoor(i) Ykoor(i) Zkoor(i)] = ecef2enu(PZ90(i,1),PZ90(i,2),PZ90(i,3),N,E,High,wgs84Ellipsoid,'radians');
     if Zkoor(i) > 0
     r(i) = sqrt(Xkoor(i)^2 + Ykoor(i)^2 + Zkoor(i)^2);
     teta(i) = acos(Zkoor(i)/r(i));
     if Xkoor(i) > 0
     phi(i) = -atan(Ykoor(i)/Xkoor(i))+pi/2;
     elseif (Xkoor(i)<0)&&(Ykoor(i)>0)
     phi(i) = -atan(Ykoor(i)/Xkoor(i))+3*pi/2;
     elseif (Xkoor(i)<0)&&(Ykoor(i)<0)
     phi(i) = -atan(Ykoor(i)/Xkoor(i))-pi/2;
     end
    else teta(i) = NaN;
     r(i) = NaN;
     phi(i) = NaN;
    end
end
%Построение графиков
figure(1)
[Xx,Yy,Zz]=sphere(30);
surf(Rz*Xx,Rz*Yy,Rz*Zz)
hold on
grid on
plot3(F1(:,1), F1(:,2), F1(:,3), 'b')
xlabel('Ось Х, м')
ylabel('Ось Y, м')
zlabel('Ось Z, м')
surf(Rz*Xx,Rz*Yy,Rz*Zz)
grid on
title({'Траектория движения спутника ГЛОНАСС №3'})
plot3(PZ90(:,1),PZ90(:,2),PZ90(:,3),'r')
xlabel('Ось Х, м')
ylabel('Ось Y, м')
zlabel('Ось Z, м')
hold off
figure (2)
pax = polaraxes;
polarplot(pax,phi,teta*180/pi,'r')
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
title('SkyView спутника ГЛОНАСС №3')
th = hours(t1./3600-3); 
figure(3);
grid on
hold on
plot(th,(-teta*180/pi+90),'DurationTickFormat','hh:mm:ss')
title('Угол места')
xlabel('Время в МДВ')
ylabel('Угол места спутника ГЛОНАСС №3, град')
