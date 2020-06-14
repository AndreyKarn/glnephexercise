% ���������, ���������� �� ������ �����
% 22 20  2 10 13 45  0.0 -.134212896228E-04 -.181898940355E-11  .495000000000E+05
%      .311761962891E+04 -.179497814178E+01 -.186264514923E-08  .000000000000E+00
%      .158781806641E+05 -.202221393585E+01 -.186264514923E-08 -.300000000000E+01
%      .196852387695E+05  .192774677277E+01 -.186264514923E-08  .000000000000E+00

clear all;
close all;

tic;
% ���������
dt = 1.0E-01;               % ��� ��������� �������� �
w_e = 7.292115E-05;      % ������� �������� �������� �����
mu = 3.986004418E+14;   % ��������� ��������������� ���� �����
R_z = 6371;              % ������ �����, ��
% ���������
T_Omega = .495000000000E+05 + 18 + 3*3600;  % ����� ������� ��������
% ����������
X = .311761962891E+07;                     
Y = .158781806641E+08;
Z = .196852387695E+08;
% ��������
Vx = -.179497814178E+04;
Vy = -.202221393585E+04;
Vz = .192774677277E+04;
% ������������ ��������� �� ������ �������� ��� (���� � ������) 
Ax = -.186264514923E-05;
Ay = -.186264514923E-05;
Az = -.186264514923E-05;
LL_potr = [55.756735, 37.703177 200];           % ���������� �����������
num = fix(12*3600/dt);                          % ���������� �������� �� 12 ������� �������� �������
num_eph = fix((T_Omega-12*3600-3*3600)/dt);     % ����� ������� � ������������ �������
t_12h = 12*60*60+3*60*60 + dt.*(1:1:num);          % ������ �������� �������
coordPZ90 = zeros(num,3);                        % ������� ��������� � ������� �����������
coordECI = zeros(num,3);
t_G0 = (9*3600+18*60+10.5009+3*3600);           %9:18:10.5009; �������� �������� ����� �� ������������� ������� ������� ����
t_G = t_G0 + w_e*(t_12h(num_eph)- 3*3600);

% �������� ��������� �� ��-90 � ECI
Xa = X*cos(t_G) - Y*sin(t_G); 
Ya = X*sin(t_G) + Y*cos(t_G);
Za = Z;
Vxa = Vx*cos(t_G) - Vy*sin(t_G) - w_e*Ya;
Vya = Vx*sin(t_G) + Vy*cos(t_G) + w_e*Xa;
Vza = Vz;
Yn = [Xa Ya Za Vxa Vya Vza];


% �������������� ������� �����-�����
% ������ ��������� �� ������ ��������� �� ������� ������� ��������
[t, Yn1] = ode45('proizv', T_Omega:-dt:t_12h(1), Yn);
coordECI(1:num_eph,:) = Yn1(end:-1:1,1:3);
% ������ ��������� �� ������� ������� �������� �� ����� ������������
% ���������
[t, Yn1] = ode45('proizv', T_Omega:dt:t_12h(num), Yn);
coordECI(num_eph:end,:) = Yn1(1:end,1:3);
% �������� ���������� ��������� �� ECI � ��-90
for i = 1:num
    t_G = t_G0 + w_e*(t_12h(i)- 3*3600);
    coordPZ90(i,1) = coordECI(i,1)*cos(t_G) + coordECI(i,2)*sin(t_G);
    coordPZ90(i,2) = -coordECI(i,1)*sin(t_G) + coordECI(i,2)*cos(t_G);
    coordPZ90(i,3) = coordECI(i,3);
end
% �������� ��������� �� ��-90 � WGS84 (��� ��������� SkyView)
ppb = 1e-9;
mas = 1e-3/206264.8; % [���]
M_WGS84 = [-3*ppb -353*mas -4*mas;
            353*mas -3*ppb 19*mas;
            4*mas -19*mas -3*ppb];
       
koord_WGS84 = coordPZ90.'; % ������� � �������-�������
for i = 1:length(koord_WGS84(1,:))
    koord_WGS84(:,i) =  koord_WGS84(:,i) + M_WGS84 * koord_WGS84(:,i) + [0.07; -0; -0.77];
end
koord_WGS84 = koord_WGS84.'; % ������� � �������-������

% �������� ��������� �� WGS84 � SkyView
X = zeros(num,1);
Y = zeros(num,1);
Z = zeros(num,1);
r = zeros(num,1);
teta = zeros(num,1);
phi = zeros(num,1);
for i = 1:length(koord_WGS84(:,1))
    [X(i), Y(i), Z(i)] = ecef2enu(koord_WGS84(i,1),koord_WGS84(i,2),koord_WGS84(i,3),LL_potr(1),LL_potr(2),LL_potr(3),wgs84Ellipsoid);
    if Z(i) > 0
        r(i) = sqrt(X(i)^2 + Y(i)^2 + Z(i)^2);
        teta(i) = acos(Z(i)/r(i));
        if X(i) > 0
            phi(i) = -atan(Y(i)/X(i))+pi/2;
        elseif (X(i)<0)&&(Y(i)>0)
            phi(i) = -atan(Y(i)/X(i))+3*pi/2;
        elseif (X(i)<0)&&(Y(i)<0)
            phi(i) = -atan(Y(i)/X(i))-pi/2;
        end
    else
        teta(i) = NaN;
        r(i) = NaN;
        phi(i) = NaN;
    end
end

% ���������� ��������

% ������ �����, ����������� �����
[Xsf, Ysf, Zsf] = sphere(25);
% ����������, ��� ����������� ���� �������
alfa1 = pi/180.*(1:359)';
beta1 = 85.*ones(359,1);

% ������ ��������� � ��-90
plot3(coordPZ90(:,1)/1000, coordPZ90(:,2)/1000, coordPZ90(:,3)/1000);
hold on
grid on
title('��������� �������� � �� ��-90');
xlabel('OX, ��');
ylabel('OY, ��');
zlabel('OZ, ��');
axis('square');
axis('equal');
surf(Xsf*R_z, Ysf*R_z, Zsf*R_z);
hold off

% ������ ��������� � ECI
figure;
plot3(coordECI(:,1)/1000, coordECI(:,2)/1000, coordECI(:,3)/1000);
hold on
grid on
title('��������� �������� � �� ECI');
xlabel('OX, ��');
ylabel('OY, ��');
zlabel('OZ, ��');
axis('square');
axis('equal');
surf(Xsf*R_z, Ysf*R_z, Zsf*R_z);
hold off

% SkyView �������� (����-�����) ������������ ������� �
figure;
ax = polaraxes;
hold on
polarplot(ax,phi,teta*180/pi,'r')
polarplot(ax,alfa1,beta1,'b')
hold off
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';
title('SkyView ������� �22')

t_12h(411564)
toc;