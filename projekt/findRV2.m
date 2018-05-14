function f_b = findRV2(e, t, p, presure, t_index ,rand)
%-------------------------------------------------------------
%  Ber?knar f_b om elementet p? index t_index i t ?r ett randelement.
%  Der?knas integralen
%  
%
% INPUT:  e,t,p,a_air, T_inf
%         rand = 3 f?r top och =1 f?r botten
%         t_index   aktuellt index i matrisen t
%
% OUTPUT: RV :      Matix 3 x 1
%-------------------------------------------------------------
point_1 = t(1,t_index);
point_2 = t(2,t_index);
point_3 = t(3,t_index);
f_b = zeros(6,1);
%det vi kallar Ne ?r inte Ne egentligen!
for i = 1:length(e)
   if (e(5,i) == rand)
       [true ,Ne] = isBoundary2(point_1, point_2, point_3,e(1,i),e(2,i));
       if(true == 1)
           p1 = p(:,e(1,i));
           p2 = p(:,e(2,i));
           Length=pdist([p1 p2]);
           s1 = sum(sum(Ne));
           
           if (rand == 3)
                ortpresure = presure * [p2(2)-p1(2); -(p2(1)-p1(1))]/Length; 
           end
           if (rand == 1)
                ortpresure = presure * [-(p2(2)-p1(2)); (p2(1)-p1(1))]/Length; 
           end
           f_b = Length * Ne * ortpresure/s1;    
           %Kce = Ne*Ne'*Length*a_air/s2;   
           
       end
    end

end