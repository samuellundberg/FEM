%% geometrin
clear all
clf
mesh = load('mesh1.mat');    %Laddar den sparade genomerin av en rektangel
e = mesh.e;                  %Randen        (edges)
p = mesh.p;                  %Noderna       (points)
t = mesh.t;                  %Elementen     (triangels)
[m,n] = size(p);

%fixa randvilkor
bc = [];
for r = 1:n
    if p(2,r) - 1.5 == 0       %fr?n labben   
        bc = [bc;r, 0];
    end
end


%Flyttar alla noder med avst?ndet -6.5*sin(pi*x/15) f?r att f? geometrin av
%den ?nskade sinusv?gen sinusv?g
for i = 1:n
    p(2,i)= p(2,i)+(-6.5*sin(pi/15*p(1,i)));
end 

%Plockar fram n?dv?ndiga parametrar som vi kommer anv?nda 

nelm=length(t(1,:));                        %Antalet element
edof(:,1)=1:nelm;                           %L?gger in elementen i Edof p? rad 1
edof(:,2:4)=t(1:3,:)';                      %L?gger in elementens noder p? rad 2,3,4
coord=p';                                  %D?per om och transponerar punkterna
ndof=max(max(t(1:3,:)));                    %Ber?knar antalet frihetsgrader
[ex,ey]=coordxtr(edof,coord,(1:ndof)',3);   %Tar fram x- och ykoordinater f?r elementen

%Plottar genometrin
%eldraw2(ex,ey,[1,2,0]);
eldraw2(ex,ey,[1 2 0])
%% Parametrar som kommer anv?ndas under projektet
ep = 1;                 %Tjocklek av geometrin
eq = 0;                 %V?rmef?rs?rjning per volymenhet

T_0=25;                 %Intiella temperaturen
T_inf = 25;             %Temperaturen "l?ngt bort"
T_ref = 50;             %Temperaturen efter tiden t_p = uppstartningstiden

t_p = 10;               %Uppstartningstiden

%Funktion som beskriver temperaturen l?ngs rand efter tiden t
% f_t = @(tt) (tt/t_p).^3.*(6.*(tt/t_p).^2-15.*tt/t_p+10);
% T_p = @(tt) (T_ref-T_0)*f_t(tt)+T_0;


a_air = 1000;           %Konvektionskoefficient 
%q_n = a_air*(t-T_inf); %Newtons konvektionslag, beskriver v?rmefl?de
alpha = 1.6e-5;         %Termisk expansionskonstant
k = 16.2;               %Konduktivitetskonstant
D = k*eye(2);           %Konduktivitetsmatrisen
rho = 7990;             %Densitet f?r v?rt material [Kgm^?3]



%% Station?ra l?sning
%ser okej ut men fel rv i ovansidan, fixas med bc
K = zeros(ndof);
f = zeros(ndof,1);
Kce = zeros(3);
a = T_0*ones(ndof,1);

for i = 1:nelm
   [Ke, fe] = flw2te(ex(i,:),ey(i,:),ep,D,eq);  
   
   %Ce=plantml(ex(i,:),ey(i,:),rho);  Skippar denna i station?ra l?sningen 
   %vi ska inte g?ra detta f?r toppen, ?ndast botten, vi ska g?ra f?r
   %toppen
   [f_b_bot, Kce_bot] = findRV(e,t,p,a_air,i,1,T_ref);         %Tar fram randvillkor om elementet vi ?r p? ?r ett randelement p? den "undre" randen
   Ke = Ke + Kce_bot;
   fe = fe + f_b_bot;
   
   [K, f] = assem(edof(i, :), K, Ke, f, fe);
end

a = a + solveq(K,f, bc);

ed=extract(edof,a);
fill(ex',ey',ed');
colorbar;
%% Numdiff
M = 100;
tend = 10;
dt = tend/10;
tt = linspace(0,tend, M+1);

for j = 1:M
    
    K = zeros(ndof);
    f = zeros(ndof,1);
    Kce = zeros(3);
    a = T_0*ones(ndof,1);
    C = zeros(ndof);
    
    for i = 1:nelm
        [Ke, fe] = flw2te(ex(i,:),ey(i,:),ep,D,eq);    
        Ce=plantml(ex(i,:),ey(i,:),rho);                                  %Skapar Ce-matrisen fr?n kurshemsidan
        [f_b_bot, Kce_bot] = findRV(e,t,p,a_air,i,1,T_p(tt(j)+dt, T_ref, T_0,t_p));         %Tar fram randvillkor om elementet vi ?r p? ?r ett randelement p? den "undre" randen
        [f_b_top, Kce_top] = findRV(e,t,p,a_air,i,3,T_inf);         %Tar fram randvillkor om elementet vi ?r p? ?r ett randelement p? den "?vre" randen
        Ke = Ke + Kce_top + Kce_bot;
        fe = fe + f_b_bot + f_b_top;
        C = assem(edof(i,:),C,Ce);
        [K, f] = assem(edof(i, :), K, Ke, f, fe);
    end

    a = a + solveq(K,f);
    a = anew(a,f,C,K,dt);
    
    ed=extract(edof,a);
    
    if(tt(j) == 2.5 ||tt(j) == 5 ||tt(j) == 7.5)
        figure
        fill(ex',ey',ed');
        colorbar;
    end
end