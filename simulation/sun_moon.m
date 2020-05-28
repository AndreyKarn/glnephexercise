function [ sun_moon ] = sun_moon( T, xa,ya,za,vxa,vya,vza )
pi=3.14159265359;
J02=1082625.75*10^-9;  %зональный гармонический коэффициент второй степени
ae=6378136;         %больша€ (экваториальна€) полуось общеземного эллипсоида
GM=398600441.8*10^6;    %Ц геоцентрическа€ константа гравитационного пол€ «емли 
J20=1082625.75*10^-9;    %зональный гармонический коэффициент второй степени, характеризующий пол€рное сжатие «емли
we=0.7292115*10^-4;        %earth's rotation rate

q_m=2.3555557435+8328.6914257190*T+ 0.0001545547*T^2;%средн€€ аномали€ Ћуны, рад
OM_m=2.1824391966-33.7570459536*T+ 0.0000362262*T^2;%средн€€ долгота восход€щего узла Ћуны
G_sh=1.4547885346+71.0176852437*T-0.0001801481*T^2;%средн€€ долгота периге€ орбиты Ћуны, рад
eps=0.4090926006-0.0002270711*T;    %средний наклон эклиптики к экватору


%sun
a_s=1.49598*10^8;%больша€ полуось Ђорбитыї —олнца
e_s=0.016719;
q_s=6.2400601269+628.3019551714*T-0.0000026820*T^2;
OM_s=-7.6281824375+0.0300101976*T+0.0000079741*T^2;


%Em находим итерационным методом, за начальное приближение берем q_m
Es=q_s; % нач приближение 
Es1=q_s-1; % нач приближение
                                                    %%ѕ–ќ¬≈–»“№!!!!!!!
while ((Es-Es1)>10e-8)
   Es1=Es;
    Es=q_s+e_s*sin(Es1);
end
TET_s=asin((sqrt(1-e_s^2)*sin(Es))/(1-e_s*cos(Es)));

eps_s=cos(TET_s)*cos(OM_s)-sin(TET_s)*sin(OM_s);
n_s=(sin(TET_s)*cos(OM_s)+cos(TET_s)*sin(OM_s))*cos(eps);
J_s=(sin(TET_s)*cos(OM_s)+cos(TET_s)*sin(OM_s))*sin(eps);
r_s=a_s*(1-e_s*cos(Es));


%moon
am= 3.84385243*10^5;%больша€ полуось орбиты Ћуны
e_m=0.054900489;%ecentrissity
i_m=0.0898041080;%rad, среднее наклонение орбиты Ћуны к плоскости эклиптики



EPS=1-(cos(OM_m))^2*OM_m*(1-cos(i_m));
nu=sin(OM_m)*sin(i_m);
J=cos(OM_m)*sin(i_m);

eps_11=sin(OM_m)*cos(OM_m)*(1-cos(i_m));
eps_12=1-(1-cos(i_m))*(sin(OM_m))^2;
n_11=EPS*cos(eps)-J*sin(eps);
n_12=eps_11*cos(eps)+nu*sin(eps);
J_11=EPS*sin(eps)+J*cos(eps);
J_12=eps_11*sin(eps)-nu*cos(eps);

%Em находим итерационным методом, за начальное приближение берем q_m
Em=q_m; % нач приближение 
Em1=q_m-1; % нач приближение
                                                    %%ѕ–ќ¬≈–»“№!!!!!!!
while ((Em-Em1)>10e-8)
   Em1=Em;
    Em=q_m+e_m*sin(Em1);
end

TET_m=asin((sqrt(1-e_m^2)*sin(Em))/(1-e_m*cos(Em)));
G_m=4902.799*10^9;%константа гравитационного пол€ Ћуны
eps_m=(sin(TET_m)*cos(G_sh)+cos(TET_m)*sin(G_sh))*eps_11+(cos(TET_m)*cos(G_sh)-sin(TET_m)*sin(G_sh))*eps_12;
n_m=(sin(TET_m)*cos(G_sh)+cos(TET_m)*sin(G_sh))*n_11+(cos(TET_m)*cos(G_sh)-sin(TET_m)*sin(G_sh))*n_12;
J_m=(sin(TET_m)*cos(G_sh)+cos(TET_m)*sin(G_sh))*J_11+(cos(TET_m)*cos(G_sh)-sin(TET_m)*sin(G_sh))*J_12;
r_m=am*(1-e_m*cos(Em));
%че за жесть
xrat_m=xa/r_m;
yrat_m=ya/r_m;
zrat_m=za/r_m;
Grat_m=G_m/r_m^2;
det_m=((eps_m-xrat_m)^2+(n_m-yrat_m)^2+(J_m-zrat_m)^2)^(3/2);

Jxam=Grat_m*(((eps_m-xrat_m)/det_m)-eps_m);
Jyam=Grat_m*(((n_m-yrat_m)/det_m)-n_m);
Jzam=Grat_m*(((J_m-zrat_m)/det_m)-J_m);

G_s=13271244.0*10^13;   %константа гравитационного пол€ —олнца
xrat_s=xa/r_s;
yrat_s=ya/r_s;
zrat_s=za/r_s;
Grat_s=G_s/r_s^2;
det_s=((eps_s-xrat_s)^2+(n_s-yrat_s)^2+(J_s-zrat_s)^2)^(3/2);

Jxas=Grat_s*(((eps_s-xrat_s)/det_s)-eps_s);
Jyas=Grat_s*(((n_s-yrat_s)/det_s)-n_s);
Jzas=Grat_s*(((J_s-zrat_s)/det_s)-J_s);

[sun_moon]=[Jxam;Jyam;Jzam;Jxas;Jyas;Jzas];
end

