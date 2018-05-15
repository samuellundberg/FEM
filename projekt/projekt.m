%% geometri & konstanter
mesh = load('mesh1.mat');    %Laddar den sparade genomerin av en rektangel
e = mesh.e;                  %Randen        (edges)
p = mesh.p;                  %Noderna       (points)
t = mesh.t;                  %Elementen     (triangels)
[m,n] = size(p);

%allokerar vektorer f?r randvilkor
bc_lower = [];
bc_upper = [];
bc_left = [];
bc_right = [];
for l = 1:n
    if p(2,l) == -1.5         
        bc_lower = [bc_lower; l, 0];
    end
    if p(2,l) == 1.5
        bc_upper = [bc_upper; l, 0];
    end
    if p(1,l) == -15.0
        bc_left = [bc_left; l, 0];
    end
    if p(1,l) == 15.0
        bc_right = [bc_right; l, 0];
    end
end

%Boj balken
for i = 1:n
    p(2,i)= p(2,i)+(-6.5*sin(pi/15*p(1,i)));
end 
p = p/1000;                                 %G�r om till meter
%Plockar fram n?dv?ndiga parametrar som vi kommer anv?nda 
nelm=length(t(1,:));                        %Antalet element
edof(:,1)=1:nelm;                           %L?gger in elementen i Edof p? rad 1
edof(:,2:4)=t(1:3,:)';                      %L?gger in elementens noder p? rad 2,3,4
coord=p';                                   %D?per om och transponerar punkterna
ndof=length(p);                             %Ber?knar antalet frihetsgrader
[ex,ey]=coordxtr(edof,coord,(1:ndof)',3);   %Tar fram x- och ykoordinater f?r elementen

%Plottar genometrin
eldraw2(ex,ey,[1 2 0],edof(:,1));
% Parametrar som kommer anv?ndas under projektet
ep = 1;                  %Tjocklek av geometrin
eq = 0;                  %V?rmef?rs?rjning per volymenhet
%begynnelsev?rden och randvilkor
T_0 = 25;                %Intiella temperaturen, samma som Tinf
T_inf = 25;
T_ref = 50;              %Temperaturen efter tiden t_p = uppstartningstiden
t_p = 10;                %Uppstartningstiden
a_air = 1000;            %Konvektionskoefficient
p_air = 10^10;
p_ref = 5*10^10;
%matrialkonstanter
%q_n = a_air*(t-T_inf);  %Newtons konvektionslag, beskriver v?rmefl?de
alpha = 1.6*10^(-5);     %Termisk expansionskonstant
k = 16.2;                %Konduktivitetskonstant
D = k*eye(2);            %Konduktivitetsmatrisen
rho = 7990;              %Densitet f?r v?rt material [Kgm^?3]
E = 1.93*10^11;          %youngs modulus
v = 0.25;                %poissons ratio
cp = 500;                %specifik v?rmekapicitet
%% Stationary solution
K = zeros(ndof);
fb = zeros(ndof,1);
bc_lower(:,2) = T_ref;

for i = 1:nelm
    Ke = flw2te(ex(i,:),ey(i,:),ep,D);
    [fb_top, Kce_top] = findRV(e, t, p, a_air, i, 3, T_inf);
    Ke = Ke + Kce_top;
    [K, fb] = assem(edof(i, :), K, Ke, fb, fb_top);
end

f = fb;
figure
statT = solveq(K, f, bc_lower);     
ed=extract(edof,statT);
fill(ex',ey',ed');
colorbar;
%% Numdiff
M = 16;
tend = 10;
dt = tend/M;
tt = linspace(0,tend, M+1);
a = T_0*ones(ndof,1);
for j = 1:M
    K = zeros(ndof);
    C = zeros(ndof);
    fb = zeros(ndof,1);

    cur_temp = T_p(tt(j)+dt, T_ref, T_0,t_p);
    bc_lower(:,2) = cur_temp;
    
    for i = 1:nelm
        Ke = flw2te(ex(i,:), ey(i,:), ep, D);       
        Ce = plantml(ex(i,:), ey(i,:), cp*rho);    
        C = assem(edof(i,:),C,Ce);
        [fb_top, Kce_top] = findRV(e, t, p, a_air, i, 3, T_0);
        Ke = Ke + Kce_top;
        [K, fb] = assem(edof(i, :), K, Ke, fb, fb_top);
    end
    
    f = fb;
    a = solveq(C+K*dt, C*a + f*dt, bc_lower);
    if(tt(j) == 2.5 ||tt(j) == 5 ||tt(j) == 7.5) 
        figure
        ed=extract(edof,a);
        fill(ex',ey',ed');
        colorbar;
    end
end

%% Stress

ptype = 2;                     %plain strain
ep = [ptype, 1];               %ep = [ptype t], ptype: 1=plain stress, 2=plain strain
D_hooke = hooke(ptype, E, v);
eq = [0, 0];             
%dT = sum(statT)/ndof-T_0;     %medeltemperaturf?r?ndringen, bra formel enligt Matias
%epsilon0 = alpha*dT*[1 1 0]'; %formel kap12, plain strain
bc_left(:,2) = 0;             %sitter fast i kanterna
bc_right(:,2) = 0;
bc = [bc_left; bc_right];
K = zeros(2*ndof);            %stiffnes,      plante.m
fb = zeros(2*ndof,1);         %kraft p? rand, r?kna ut sj?lv..
f0 =  zeros(2*ndof,1);           %termisk last,  plantf.m + snake

ex_wide = [ex(:,1) ex(:,1) ex(:,2) ex(:,2) ex(:,3) ex(:,3)];
ey_wide = [ey(:,1) ey(:,1) ey(:,2) ey(:,2) ey(:,3) ey(:,3)];
%ep = [ptype t ], ptype: analysis type, t: thickness, kolumn 1 ?r ny..
%det vi m?ste g?ra ?r att skapa en ny edof, vi ska g? fr?n 1 till 2 
%frihetsgrader per element, det ?r lite oklart hur vi ska hitta de nya
%v?rdena, typ ta ggr2 eller n?got..

%Vi beh?ver expandera edof, d? varje nod nu har tv? d.o.f.

c1 = [edof(:,2)*2 - 1, edof(:,2)*2];
c2 = [edof(:,3)*2 - 1, edof(:,3)*2];
c3 = [edof(:,4)*2 - 1, edof(:,4)*2];
edof_wide = [edof(:,1), c1, c2, c3];
for i = 1:nelm
    Ke = plante(ex(i,:), ey(i,:), ep, D_hooke);
    dT = (statT(t(1,i))+statT(t(2,i))+statT(t(3,i)))/3 - T_0;
    epsilon0 = alpha*dT*[1 1 0]';
    fe0 = plantf(ex(i,:), ey(i,:), ep, (D_hooke([1 2 4], [1 2 4])*epsilon0)');
    fb_top = findRV2(e, t, p, p_air, i ,3);
    fb_bot = findRV2(e, t, p, p_ref, i ,1);
    fe0 = fe0 + fb_top + fb_bot;
    [K, f0] = assem(edof_wide(i, :), K, Ke, f0, fe0);
end
f=f0;

u = solveq(K, f, bc);  
ed=extract(edof_wide,u);

% VON MISES
Seff_el = zeros(nelm, 1);
Seff_nod = zeros(ndof, 1);
for i = 1:nelm
    [es, et] = plants(ex(i,:), ey(i,:), ep, D_hooke([1 2 4], [1 2 4]), ed(i,:));
    %es1(i,:) = es;
    sigma_zz=et(1)*D(2,1); %Sida 315 i kursboken, D(2,1) tar ut E*v/((1+v)(1-2v)) och et(1) = epsilon_xx
    Mises = sqrt(es(1,1)^2 + es(1,2)^2 + sigma_zz^2 - es(1,1)*es(1,2) - es(1,1)*sigma_zz - es(1,2)*sigma_zz + 3*es(1,3)^2);
    Seff_el(i,1) = Mises;
end

for i=1:ndof
    [c0,c1]=find(edof(:,2:4)==i);
    Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1);
end

figure(1)
ed1 = extract(edof, Seff_nod);
fill(ex',ey',ed1')
colorbar;

figure(2)
eldraw2(ex,ey, [1 2 1]);
eldisp2(ex,ey, ed, [1 4 1], 100);