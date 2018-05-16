function [fb, Kce] = findRV(e, t, p, a_air, t_index ,rand, T_inf)
%-------------------------------------------------------------
%  Ber?knar f_b om elementet p? index t_index i t ?r ett randelement.
%  Der?knas integralen
%  a_air*T_inf*int(N')
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
fb = zeros(3,1);
Kce = zeros(3);
%det vi kallar Ne ?r inte Ne egentligen!
for i = 1:length(e)
   if (e(5,i) == 3)
        [a, Ne, Ne2] = isBoundary(point_1, point_2, point_3, e(1,i), e(2,i));
        if (a == 1)
            p1 = p(:,e(1,i));
            p2 = p(:,e(2,i));
            Length=sqrt((p1(1,1)-p2(1,1)).^2+(p1(2,1)-p2(2,1)).^2);           
%             s1 = sum(Ne);
%             s2 = sum(sum(Ne2));
            fb = Length * a_air * T_inf * Ne /2;
            Kce = Ne2 * Length * a_air /6;
            if (isequal(p1,[0.015;0.0015]) || isequal(p2,[0.015;0.0015]))
               fb = Length * a_air * T_inf * [2;1;0] /2;
            end
            if (isequal(p1,[-0.015;0.0015]) || isequal(p2,[-0.015;0.0015]))
               fb = Length * a_air * T_inf * [2;1;0] /2; 
            end

        end
   end
end

end