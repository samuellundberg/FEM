function [a, Ne] = isBoundary2(point_1, point_2, point_3, Bp1, Bp2)
%Jämför om någon av punkterna p1 p2 p3 överrensstämmer med 
%Randpunkterna Bp1 och Bp2
if(point_1 == Bp1 && point_2 == Bp2)
    Ne = [1 0 1 0 0 0;
          0 1 0 1 0 0]';
    a = 1;
elseif(point_2 == Bp1 && point_1 == Bp2)
    Ne = [1 0 1 0 0 0;
          0 1 0 1 0 0]';
    a = 1;
elseif(point_1 == Bp1 && point_3 == Bp2)
    Ne = [1 0 0 0 1 0;
          0 1 0 0 0 1]';
    a = 1;
elseif(point_3 == Bp1 && point_1 == Bp2)
    Ne = [1 0 0 0 1 0;
          0 1 0 0 0 1]';
    a = 1;
elseif(point_2 == Bp1 && point_3 == Bp2)
    Ne = [0 0 1 0 1 0;
          0 0 0 1 0 1]';
    a = 1;
elseif(point_3 == Bp1 && point_2 == Bp2)
    Ne = [0 0 1 0 1 0;
          0 0 0 1 0 1]';
    a = 1;
else
    Ne = [0 0 0 0 0 0;
          0 0 0 0 0 0]';
    a = 0;
end



end