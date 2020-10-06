function dF = diffs(t, F)
%константы
GM = 398600441.8e6;
J02 = -1082.63*10^-6;
ae = 6378136;
%расчет
Xa=F(1);
Ya=F(2);
Za=F(3);
r=sqrt(Xa^2 + Ya^2 + Za^2);
GMn = GM/r^2;
x1 = Xa/r;
y1 = Ya/r;
z1 = Za/r;
RO = ae/r;
dF = F(:);
dF(1) = F(4);
dF(2) = F(5);
dF(3) = F(6); 
dF(4) = -GMn*x1 + 1.5*J02*GMn*x1*(RO^2)*(1 - 5*z1^2);
dF(5) = -GMn*y1 + 1.5*J02*GMn*y1*(RO^2)*(1 - 5*z1^2);
dF(6) = -GMn*z1 + 1.5*J02*GMn*z1*(RO^2)*(3 - 5*z1^2);
end

