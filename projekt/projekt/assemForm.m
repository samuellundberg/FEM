function fbvec = assemForm(t, p, index, fbvec,rand)
%-------------------------------------------------------------
%  Calculates the form fucktion for the bondaryelement at i 
%  puts it at fbvec_top(:,index)
%  
% INPUT:  t,p, index, fbvec_top
%         rand = 3 -> top
%         rand = 1 -> bot
% OUTPUT: RV :      Matix 3 x nelm
%-------------------------------------------------------------

elenod = t([1 2 3],index);
if rand==3
    pos = 1.5;
end
if rand==1
    pos = -1.5;
end
if(p(2,elenod(1))== pos && p(2,elenod(2)) == pos)
    form = [1; 1; 0];
    fbvec(:,index) = form; 
elseif(p(2,elenod(1))== pos && p(2,elenod(3)) == pos)
    form = [1; 0; 1];
    fbvec(:,index) = form;
elseif(p(2,elenod(2))== pos && p(2,elenod(3)) == pos)
    form = [0; 1; 1];
    fbvec(:,index) = form;
end
end