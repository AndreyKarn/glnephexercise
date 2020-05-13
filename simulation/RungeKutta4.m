function [t, Y] = RungeKutta4( tn, h, Tlim, Y0 )

t = tn:h:Tlim;
t = t.';

Y(1,:) = Y0;
    
% KnextstepC = load('../test_data_out.txt');

for k = 2:length(t)
    
    K1 = diffs(tn, Y(k-1,:));
    
    Y2 = Y(k-1,:) + h*K1.'/2;
    K2 = diffs(tn + h/2, Y2);
    
    Y3 = Y(k-1,:) + h*K2.'/2;
    K3 = diffs(tn + h/2, Y3);
    
    Y4 = Y(k-1,:) + h*K3.';
    K4 = diffs(tn + h, Y4);
    
    Knextstep = h/6 * (K1 + 2*K2 + 2*K3 + K4);
    Y(k,:) = Y(k-1,:) + Knextstep.';
    
%     deltaK11(k) = KnextstepC(k-1,1) - Knextstep(1);
%     deltaK12(k) = KnextstepC(k-1,2) - Knextstep(2);
%     deltaK13(k) = KnextstepC(k-1,3) - Knextstep(3);
%     deltaK14(k) = KnextstepC(k-1,4) - Knextstep(4);
%     deltaK15(k) = KnextstepC(k-1,5) - Knextstep(5);
%     deltaK16(k) = KnextstepC(k-1,6) - Knextstep(6);
end

% figure
% hold on
% plot(t, deltaK11, t, deltaK12, t, deltaK13, t, deltaK14, t, deltaK15, t, deltaK16)
% legend('deltaK11','deltaK12','deltaK13','deltaK14','deltaK15','deltaK16')

end