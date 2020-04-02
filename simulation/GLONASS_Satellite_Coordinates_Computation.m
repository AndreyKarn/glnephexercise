clear all; close all; tic; clc;
format long;

%% Заданные параметры
% Ephemerides reference epoch:
% 2020/02/25 13:45:18
Te = 18 + 45*60 + 13*60*60; % 49518 s

% Coordinate at Te in PZ-90:
X = -8444572.27; % m
Y = -8664957.52; % m
Z = 22466454.10; % m

% Velocity component at Te in PZ-90:
VX = 2983.60348; % m/s
VY = -743.76965; % m/s
VZ =  832.83615; % m/s

% Moon and sun acceleration at Te:
AX = -0.0000028; % m/s2
AY =  0.0000019; % m/s2
AZ = -0.0000019; % m/s2

% SV clock offset:
Tau = -44762.2; % ns
% SV relative frequency offset:
Gamma = 0.0009; % ns/s

%% 1. Coordinates transformation to an inertial reference frame:
% earth's rotation rate:
Omega_E = 0.7292115e-4; % rad/s

%Xa_i = nan(size(900));

Toe = 12*60*60;
Tof = 24*60*60;

for k = Toe:Tof
    
    % the sidereal time in Greenwich at midnight GMT of a date at which the
    % epoch Te is specified:
    Theta_G0 = k;
    % the sidereal time at epoch Te, to which are referred the initial
    % conditions, in Greenwich meridian:
    Theta_Ge = Theta_G0 + Omega_E * (Te - 3 * 60 * 60);
    
    if k == Toe
        % Position:
        Xa = X * cos(Theta_Ge) - Y * sin(Theta_Ge);
        Ya = X * sin(Theta_Ge) + Y * cos(Theta_Ge);
        Za = Z;
        
        % Velocity:
        VXa = VX * cos(Theta_Ge) - VY * sin(Theta_Ge) - Omega_E * Ya;
        VYa = VX * sin(Theta_Ge) + VY * cos(Theta_Ge) + Omega_E * Xa;
        VZa = VZ;
        
        % Accelerations:
        JXams = AX * cos(Theta_Ge) - AY * sin(Theta_Ge);
        JYams = AX * sin(Theta_Ge) + AY * cos(Theta_Ge);
        JZams = AZ;
    else
        % Position:
        Xa = X * cos(Theta_Ge) - Y * sin(Theta_Ge);
        Ya = X * sin(Theta_Ge) + Y * cos(Theta_Ge);
        Za = Z;
        
        % Velocity:
        VXa = VXa + dXadt;
        VYa = VYa + dYadt;
        VZa = VZa + dZadt;
        
        %Accelerations:
        JXams = JXams + dVXadt;
        JYams = JYams + dVYadt;
        JZams = JZams + dVZadt;
    end
    %% 2. Numerical integration of differential equations that describe the motion of the satellites.
    % Equatorial radius of the Earth (PZ-90):
    aE = 6378.136; % km
    
    % Gravitational constant (PZ-90):
    Mu = 398600.44; % km^3/s^2
    
    % Second zonal coefficient of spherical harmonic expression:
    C20 = -1082.63e-6; %
    
    r = sqrt(Xa^2 + Ya^2 + Za^2);
    
    Mu_norm = Mu / (r^2);
    Xa_norm = Xa / r;
    Ya_norm = Ya / r;
    Za_norm = Za / r;
    Xa_norm = Xa / r;
    Rho = aE / r;
    
    dXadt = VXa;
    dYadt = VYa;
    dZadt = VZa;
    
    dVXadt = - Mu_norm * Xa_norm ...
        + 3/2 * C20 * Mu_norm * Xa_norm * Rho^2 * (1 - 5 * Za_norm^2)...
        + JXams;
    dVYadt = - Mu_norm * Ya_norm ...
        + 3/2 * C20 * Mu_norm * Ya_norm * Rho^2 * (1 - 5 * Za_norm^2)...
        + JYams;
    dVZadt = - Mu_norm * Za_norm ...
        + 3/2 * C20 * Mu_norm * Za_norm * Rho^2 * (3 - 5 * Za_norm^2)...
        + JZams;
    
    Xm(k-Toe+1) = Xa;
    Ym(k-Toe+1) = Ya;
    Zm(k-Toe+1) = Za;
end

T = Toe:Tof;
R_Earth = 6371e3;

%% построение графиков
[X,Y,Z] = sphere(30);
surf(X*R_Earth,Y*R_Earth,Z*R_Earth)
hold on
grid on
plot3(Xm, Ym, Zm, 'b')












