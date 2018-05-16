function [Ne, Ne2] = isBoundaryH(point_1, point_2, point_3, Bp1, Bp2)
%Jämför om n?gon av punkterna p1 p2 p3 överrensstämmer med 
%Randpunkterna Bp1 och Bp2
Ne  =zeros(3,1);
Ne2= zeros(3);
if (Bp1 == point_1 || Bp2 == point_1)
    Ne(1) = 1;
elseif (Bp1 == point_2 || Bp2 == point_2)
    Ne(2) = 1;
elseif (Bp1 == point_3 || Bp2 == point_3)
    Ne(3) = 1;
end

if (isequal(Ne,[1; 1; 0]))
    Ne2 = [2 1 0; 1 2 0; 0 0 0];
elseif(isequal(Ne,[1; 0; 1]))
    Ne2 = [2 0 1; 0 0 0; 1 0 2];
elseif(isequal(Ne,[0; 1; 1]))
    Ne2 = [0 0 0; 0 2 1; 0 1 2];  
end 

end