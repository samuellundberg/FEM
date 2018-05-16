function fbvec_top = assemformTop(t, p, nelm, index, fbvec_top)
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
%det vi kallar Ne ?r inte Ne egentligen!

elenod = t([1 2 3],index);
%form = zeros(3,1);
if(p(2,elenod(1))==1.5 && p(2,elenod(2)) == 1.5)
    form = [1; 1; 0];
    fbvec_top(:,index) = form;
    %Length=sqrt((p(1,elenod(2))-p(1,elenod(1))).^2+(p(2,elenod(2))-p(2,elenod(1))).^2)  
elseif(p(2,elenod(1))==1.5 && p(2,elenod(3)) == 1.5)
    form = [1; 0; 1];
    fbvec_top(:,index) = form;
    %Length=sqrt((p(1,elenod(3))-p(1,elenod(1))).^2+(p(2,elenod(3))-p(2,elenod(1))).^2);  
elseif(p(2,elenod(2))==1.5 && p(2,elenod(3)) == 1.5)
    form = [0; 1; 1];
    fbvec_top(:,index) = form;
    %Length=sqrt((p(1,elenod(2))-p(1,elenod(3))).^2+(p(2,elenod(2))-p(2,elenod(3))).^2);  
end
end