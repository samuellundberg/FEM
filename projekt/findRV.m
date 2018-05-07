function [f_b, Kce] = findRV(e, t, p, a_air, T_inf, t_index ,rand, randvillkor)
%-------------------------------------------------------------
%  Beräknar f_b om elementet på index t_index i t är ett randelement.
%  Deräknas integralen
%  a_air*T_inf*int(N')
%
% INPUT:  e,t,p,a_air, T_inf
%         rand = 3 för top och =1 för botten
%         t_index   aktuellt index i matrisen t
%
% OUTPUT: RV :      Matix 3 x 1
%-------------------------------------------------------------
point_1 = t(1,t_index);
point_2 = t(2,t_index);
point_3 = t(3,t_index);
f_b = zeros(3,1);
Kce = zeros(3);

for i = 1:length(e)
   if (e(5,i) == rand)
       [true ,Ne] = isBoundary(point_1, point_2, point_3,e(1,i),e(2,i));
       if(true == 1)
           p1 = p(:,e(1,i));
           p2 = p(:,e(2,i));
           Length=pdist([p1 p2]);
           f_b = Length*a_air*T_inf * Ne * randvillkor;
           
           Kce = Ne*Ne'*Length*a_air;
       end
   end
end



end

