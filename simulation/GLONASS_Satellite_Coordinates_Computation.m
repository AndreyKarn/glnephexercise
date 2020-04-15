clear all; close all; tic; clc;
format long;

%% �������� ���������
% ��������� �� ��������  �����:
% 2020/02/25 13:45:18
Time_year = 2020;
Time_month = 2;
Time_day = 25;
Time_hour = 13;
Time_minutes = 45;
Time_seconds = 18;

% ���������� �� Te � ������� ��-90, [�]:
X = -8444572.27;
F = -8664957.52;
Z = 22466454.10;

% ���������� ������� �������� �� Te � ������� ��-90, [�/�]:
VX = 2983.60348;
VY = -743.76965;
VZ =  832.83615;

% ��������� �����-��������� �� Te � ������� ��-90, [�/�2]:
AX = -0.0000028;
AY =  0.0000019;
AZ = -0.0000019;

% SV ��������� ��������, [��]:
Tau = -44762.2;
% SV ������������� �������� �������, [��/�]:
Gamma = 0.0009;

%% ������ ������� ������� �������
N4 = floor((Time_year-1996)/4) + 1; % ����� �������������� ���������
NT = 365*(Time_year-1996-4*(N4-1)) + 31 + Time_day + 1; % ����� ����� � ������������� ���������
tb = Time_seconds + Time_minutes*60 + Time_hour*60*60 + 10800;

% ������ �������� ��������� ������� �� ��������
GMST = GMST_calc( N4,NT );

%% �������� ��������� � ����������� ������� �������� ������ ���� ��� � ��������� � ������ ������� ��������� ��-90
% ������� ������� �������� �������� ����� ������������ ����� ��������� �������������, [���/�]:
Omega_E = 7.2921151467e-5;

Theta_Ge = GMST + Omega_E * (tb - 3 * 60 * 60);

% ����������:
X0 = X * cos(Theta_Ge) - F * sin(Theta_Ge);
Y0 = X * sin(Theta_Ge) + F * cos(Theta_Ge);
Z0 = Z;

% ��������:
VX0 = VX * cos(Theta_Ge) - VY * sin(Theta_Ge) - Omega_E * Y0;
VY0 = VX * sin(Theta_Ge) + VY * cos(Theta_Ge) + Omega_E * X0;
VZ0 = VZ;

% ���������:
JX0ms = AX * cos(Theta_Ge) - AY * sin(Theta_Ge);
JY0ms = AX * sin(Theta_Ge) + AY * cos(Theta_Ge);
JZ0ms = AZ;
        
%% �������������� ��������� �������
Toe = 12*60*60;
Tof = 24*60*60;
Ts = 0.1;

ti = tb:Ts:Tof;

%F0 = [X0 Y0 Z0 VX0 VY0 VZ0]; % ��������� �������
F0(1,:) = [X0 Y0 Z0 JX0ms JY0ms JZ0ms]; % ��������� �������

%for k = 1:length(ti+1)
    [t, F] = ode45('diffs', ti, F0(1,:));
    %F0(k+1,:) = F(k+1,:);
%end

R_Earth = 6371e3;

%% ���������� ��������
[Xz,Yz,Zz] = sphere(30);
surf(Xz*R_Earth,Yz*R_Earth,Zz*R_Earth)
hold on
grid on
plot3(F(:,1), F(:,2), F(:,3), 'b')
xlabel('��� �, �');
ylabel('��� Y, �');
zlabel('��� Z, �');

figure
plot(t,F(:,1))



