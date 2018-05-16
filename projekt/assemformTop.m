function fbvec_top = assemForm(t, p, index, fbvec_top,rand)
%-------------------------------------------------------------
%  Tar fram formfunktion för randelementetet index och lägger
%  den på platsen fbvec_top(:,index)
%  
% INPUT:  t,p, index, fbvec_top
%         rand = 1 -> top
%         rand = 2 -> bot
% OUTPUT: RV :      Matix 3 x nelm
%-------------------------------------------------------------

elenod = t([1 2 3],index);
if rand==1
    pos = 1.5;
end
if rand==2
    pos = -1.5;
end
%form = zeros(3,1);
if(p(2,elenod(1))== pos && p(2,elenod(2)) == pos)
    form = [1; 1; 0];
    fbvec_top(:,index) = form;
    %Length=sqrt((p(1,elenod(2))-p(1,elenod(1))).^2+(p(2,elenod(2))-p(2,elenod(1))).^2)  
elseif(p(2,elenod(1))==pos && p(2,elenod(3)) == pos)
    form = [1; 0; 1];
    fbvec_top(:,index) = form;
    %Length=sqrt((p(1,elenod(3))-p(1,elenod(1))).^2+(p(2,elenod(3))-p(2,elenod(1))).^2);  
elseif(p(2,elenod(2))==pos && p(2,elenod(3)) == pos)
    form = [0; 1; 1];
    fbvec_top(:,index) = form;
    %Length=sqrt((p(1,elenod(2))-p(1,elenod(3))).^2+(p(2,elenod(2))-p(2,elenod(3))).^2);  
end
end