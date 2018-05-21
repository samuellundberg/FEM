function [Fbc] = assemFbcStrain(p, t, fbvec, index, presure, rand)
%-------------------------------------------------------------
%  Calculates the form function with 2 dof for the bondary element index for 
%  Calculates the integral:
%  int(N'*presure*t)
%
% INPUT:  e,t,p,presure
%         rand = 3 for upper bondary och =1 for lower bondary
%         index is the current element in the matrix t
%
% OUTPUT: RV :      Matix 3 x 1
%-------------------------------------------------------------
Fbc = zeros(6,1);
form = fbvec(:, index);
if(form(1)~=0 || form(2)~=0  || form(3)~=0)  
    if (form(3)==0)
        p1=[p(1,t(1, index));p(2,t(1, index))];
        p2=[p(1,t(2, index));p(2,t(2, index))];
        Length=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2); 
        
        form2 = [1 0 1 0 0 0;
                 0 1 0 1 0 0]';
    elseif(form(2) == 0)
        p1=[p(1,t(1, index));p(2,t(1, index))];
        p2=[p(1,t(3, index));p(2,t(3, index))];
        Length=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2); 
        
        form2 = [1 0 0 0 1 0;
                 0 1 0 0 0 1]';
    else 
        p1=[p(1,t(2, index));p(2,t(2, index))];
        p2=[p(1,t(3, index));p(2,t(3, index))];
        Length=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2); 
        
        form2 = [0 0 1 0 1 0;
                 0 0 0 1 0 1]';
    end
    if (rand == 3)
        norm = [p2(2)-p1(2); -(p2(1)-p1(1))]/Length;
        if (norm(2)>0)
            norm = -norm; 
        end
    end
    if (rand == 1)
        norm = [-(p2(2)-p1(2)); (p2(1)-p1(1))]/Length;
        if (norm(2)<0)
            norm = -norm; 
        end 
    end
   
    Fbc = Length * presure * form2 / 2 * norm;
end
end