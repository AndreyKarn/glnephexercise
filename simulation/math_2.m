function [res] = math_2(xa,ya,za,vxa,vya,vza,Jsm_x,Jsm_y,Jsm_z,time_start,time_final, te,T)
%надо рассмотреть два случая: когда наблюдения есть до начала старта, и
%когда они есть после старта. Во втором случае задача разбивается на 2
%части: поиск координат от старта до наблюдений и от наблюдений до финала.
dt=1;
if(te<=time_start&&te<time_final)
   t=time_start:dt:(time_final-1);
    result=nan(length(t),6);
    result(1,:)=[xa,ya,za,vxa,vya,vza]; 
    result=RungKUTT( t, result, T ); 
   
    
elseif(te>time_start&&te<time_final)
    t_bef=te:-dt:time_start; %из-за того что в этом времени матрица записывается как бы наоборот, вращаем столбцы!!!!
    t_after=te:dt:(time_final-1);
    %before
      
        result_before=nan(length(t_bef),6);
        result_before(1,:)=[xa,ya,za,vxa,vya,vza]; 
        result_before=RungKUTT( t_bef, result_before, -dt ); 
        t_before=rot90(t_bef,2);
        %result_bef=RungKUTT( t_bef,result_bef, T ); 
        %вот тут вращаем, и время тоже! потому что иначе будет [te t-1 t-2
        %... t_start], и это нормально не скоеитть! а нам надо [t_start
        %t_start+1 ... te] для этого поворачиваем каждый столбец на 180
        %градусов!
        
        
               result_before = [rot90(result_before(:,1),2) rot90(result_before(:,2),2) rot90(result_before(:,3),2) rot90(result_before(:,4),2) rot90(result_before(:,5),2) rot90(result_before(:,6),2)];
    %after
        
        result_after=nan(length(t_after),6);
        result_after(1,:)=[xa,ya,za,vxa,vya,vza]; 
        result_after=RungKUTT( t_after, result_after, dt ); 
   
       
    %соединяем
       result=[result_before;result_after];
        t=[t_before t_after];


end
%%поправка на небесные тела
tau=t-te;
tau=rot90(tau);


dx=Jsm_x*0.5*tau.^2;
dy=Jsm_y*0.5*tau.^2;
dz=Jsm_z*0.5*tau.^2;

dvx=Jsm_x*tau;
dvy=Jsm_y*tau;
dvz=Jsm_z*tau;


result(:,1)=result(:,1)+dx;
result(:,2)=result(:,2)+dy;
result(:,3)=-1*(result(:,3)+dz);     % чтоб земля на место встала

result(:,4)=result(:,4)+dvx;
result(:,5)=result(:,5)+dvy;
result(:,6)=result(:,6)+dvz;
res = [result rot90(t,3)];
end

