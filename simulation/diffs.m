function dF = diffs( t, F )
%% Учет ускорений
% Ускорения лунно-солнечные на Te в системе ПЗ-90, [м/с2]:
AX = -0.0000028;
AY =  0.0000019;
AZ = -0.0000019;

Omega_E = 7.2921151467e-5;

N4 = 7; % Номер четырехлетнего интервала
NT = 57; % Номер суток в четырехлетнем интервале

if t >= 24*60*60
    NT = NT + 1;
    t = t - 24*60*60;
end

GMST = GMST_calc( N4,NT );

Theta_Ge = GMST + Omega_E * (t - 3 * 60 * 60);

JX0ms = AX * cos(Theta_Ge) - AY * sin(Theta_Ge);
JY0ms = AX * sin(Theta_Ge) + AY * cos(Theta_Ge);
JZ0ms = AZ;

%% Расчет переменных
J02 = 1082625.75e-9; % зональный гармонический коэффициент второй степени, характеризующий полярное сжатие Земли
GM = 398600441.8e6; % геоцентрическая константа гравитационного поля Земли с учетом атмосферы, [м3/c2]
a_e = 6378136; % большая полуось общеземного эллипсоида, [м]

crdX = F(1);
crdY = F(2);
crdZ = F(3);

r = sqrt(crdX^2 + crdY^2 + crdZ^2);

GM0 = GM / (r^2);
Rho = a_e / r;
crdX0 = crdX / r;
crdY0 = crdY / r;
crdZ0 = crdZ / r;

%% Дифуры
dF = F(:);
dF(1) = F(4);
dF(2) = F(5);
dF(3) = F(6);

dF(4) = - GM0 * crdX0 - 3/2 * J02 * GM0 * crdX0 * Rho^2 * (1 - 5 * crdZ0^2) + JX0ms;
dF(5) = - GM0 * crdY0 - 3/2 * J02 * GM0 * crdY0 * Rho^2 * (1 - 5 * crdZ0^2) + JY0ms;
dF(6) = - GM0 * crdZ0 - 3/2 * J02 * GM0 * crdZ0 * Rho^2 * (3 - 5 * crdZ0^2) + JZ0ms;
end