function [t, Y] = RungeKutta4( tn, h, Tlim, Y0 )

t = tn:h:Tlim;
t = t.';

Y(1,:) = Y0;
    
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
end

end