function [RungKUTT ] = RungKUTT( t, result, dt)


for i=1:length(t)-1
       
 
    K0=F(result(i,:));
  
    
    K1=F((result(i,:)+0.5*dt.*K0));
    
    K2=F(result(i,:)+0.5*dt.*K1);
    
    K3=F(result(i,:)+dt.*K2);
    
    result(i+1,:)=result(i,:)+(dt/6)*(K0+2*K1+2*K2+K3);
%     (dt/6)*(K0+2*K1+2*K2+K3)
%     K0
%     2*K1
%     2*K2
%     K3
    
end
[RungKUTT]=result;
end

