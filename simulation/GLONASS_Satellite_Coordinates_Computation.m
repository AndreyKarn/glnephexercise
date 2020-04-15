clear all; close all; tic; clc;
format long;

%% Заданные параметры
% Эфемериды на заданную  эпоху:
% 2020/02/25 13:45:18
Time_year = 2020;
Time_month = 2;
Time_day = 25;
Time_hour = 13;
Time_minutes = 45;
Time_seconds = 18;

% Координаты на Te в системе ПЗ-90, [м]:
X = -8444572.27;
F = -8664957.52;
Z = 22466454.10;

% Компоненты вектора скорости на Te в системе ПЗ-90, [м/с]:
VX = 2983.60348;
VY = -743.76965;
VZ =  832.83615;

% Ускорения лунно-солнечные на Te в системе ПЗ-90, [м/с2]:
AX = -0.0000028;
AY =  0.0000019;
AZ = -0.0000019;

% SV временное смещение, [нс]:
Tau = -44762.2;
% SV относительное смещение частоты, [нс/с]:
Gamma = 0.0009;

%% Расчет времени формата ГЛОНАСС
N4 = floor((Time_year-1996)/4) + 1; % Номер четырехлетнего интервала
NT = 365*(Time_year-1996-4*(N4-1)) + 31 + Time_day + 1; % Номер суток в четырехлетнем интервале
tb = Time_seconds + Time_minutes*60 + Time_hour*60*60 + 10800;

% Расчет среднего звездного времени по Гринвичу
GMST = GMST_calc( N4,NT );

%% Пересчет координат и оставляющих вектора скорости центра масс НКА в связанной с Землей системе координат ПЗ-90
% средняя угловая скорость вращения Земли относительно точки весеннего равноденствия, [рад/с]:
Omega_E = 7.2921151467e-5;

Theta_Ge = GMST + Omega_E * (tb - 3 * 60 * 60);

% Координаты:
X0 = X * cos(Theta_Ge) - F * sin(Theta_Ge);
Y0 = X * sin(Theta_Ge) + F * cos(Theta_Ge);
Z0 = Z;

% Скорости:
VX0 = VX * cos(Theta_Ge) - VY * sin(Theta_Ge) - Omega_E * Y0;
VY0 = VX * sin(Theta_Ge) + VY * cos(Theta_Ge) + Omega_E * X0;
VZ0 = VZ;

% Ускорения:
JX0ms = AX * cos(Theta_Ge) - AY * sin(Theta_Ge);
JY0ms = AX * sin(Theta_Ge) + AY * cos(Theta_Ge);
JZ0ms = AZ;
        
%% Интегрирование численным методом
Toe = 12*60*60;
Tof = 24*60*60;
Ts = 0.1;

ti = tb:Ts:Tof;

%F0 = [X0 Y0 Z0 VX0 VY0 VZ0]; % Начальные условия
F0(1,:) = [X0 Y0 Z0 JX0ms JY0ms JZ0ms]; % Начальные условия

%for k = 1:length(ti+1)
    [t, F] = ode45('diffs', ti, F0(1,:));
    %F0(k+1,:) = F(k+1,:);
%end

R_Earth = 6371e3;

%% построение графиков
[Xz,Yz,Zz] = sphere(30);
surf(Xz*R_Earth,Yz*R_Earth,Zz*R_Earth)
hold on
grid on
plot3(F(:,1), F(:,2), F(:,3), 'b')
xlabel('Ось Х, м');
ylabel('Ось Y, м');
zlabel('Ось Z, м');

figure
plot(t,F(:,1))



