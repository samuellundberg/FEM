function [Fbc, Kce] = assemFbcKec(p, t, fbvec, index, a_air, T_inf)
%-------------------------------------------------------------
%  Calculates the bondary vector fbc and stiffness matrix Kc.
%  The following integrals are calculated:
%  int(N'*a_air*T_inf)
%  int(N'*N *a_air)
%
% INPUT:  e,t,p,a_air, T_inf
%         rand = 3 for upper bondary och =1 for lower bondary
%         index is the current element in the matrix t
%
% OUTPUT: RV :      Matix 3 x 1
%-------------------------------------------------------------
Fbc = zeros(3,1);
Kce = zeros(3);
form = fbvec([1 2 3], index);
if(form(1)~=0 || form(2)~=0  || form(3)~=0)  
    if (form(3)==0)
        p1=[p(1,t(1, index));p(2,t(1, index))];
        p2=[p(1,t(2, index));p(2,t(2, index))];
        Length=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2); 
        
        form2 = [2 1 0; 1 2 0; 0 0 0];
    elseif(form(2) == 0)
        p1=[p(1,t(1, index));p(2,t(1, index))];
        p2=[p(1,t(3, index));p(2,t(3, index))];
        Length=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2); 
        
        form2 = [2 0 1; 0 0 0; 1 0 2];
    else 
        p1=[p(1,t(2, index));p(2,t(2, index))];
        p2=[p(1,t(3, index));p(2,t(3, index))];
        Length=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2); 
        
        form2 = [0 0 0; 0 2 1; 0 1 2];
    end
    Fbc = Length * a_air * T_inf * form /2;
    Kce = form2 * Length * a_air /6;  
end
end