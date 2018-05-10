%% geometrin

mesh = load('mesh1.mat');    %Laddar den sparade genomerin av en rektangel
e = mesh.e;                  %Randen        (edges)
p = mesh.p;                  %Noderna       (points)
t = mesh.t;                  %Elementen     (triangels)
[m,n] = size(p);
%allokerar vektorer f?r randvilkor

% bc_upper = [];  % anv?nds inte
% for r = 1:n
%     if p(2,r) - 1.5 == 0        
%         bc_upper = [bc_upper;r, 0];
%     end
% end
bc_lower = [];
for l = 1:n
    if p(2,l) == -1.5       %fr?n labben   
        bc_lower = [bc_lower; l, 0];
    end
end
%B?j balken
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
eldraw2(ex,ey,[1 2 0]);
% Parametrar som kommer anv?ndas under projektet
ep = 1;                 %Tjocklek av geometrin
eq = 0;                 %V?rmef?rs?rjning per volymenhet
%begynnelsev?rden och randvilkor
T_0 = 25;               %Intiella temperaturen, samma som Tinf
T_ref = 50;             %Temperaturen efter tiden t_p = uppstartningstiden
t_p = 10;               %Uppstartningstiden
a_air = 1000;           %Konvektionskoefficient
p_air = 10^10;
p_ref = 5*10^10;
%matrialkonstanter
%q_n = a_air*(t-T_inf); %Newtons konvektionslag, beskriver v?rmefl?de
alpha = 1.6*10^(-5);         %Termisk expansionskonstant
k = 16.2;               %Konduktivitetskonstant
D = k*eye(2);           %Konduktivitetsmatrisen
rho = 7990;             %Densitet f?r v?rt material [Kgm^?3]
E = 1.93*10^11;          %youngs modulus
v = 0.25;                %poissons ratio
cp = 500;                %specifik v?rmekapicitet
%% Station?ra l?sning
K = zeros(ndof);
Kc = zeros(ndof);
fb = zeros(ndof,1);
bc_lower(:,2) = T_ref;
for i = 1:nelm
    Ke = flw2te(ex(i,:),ey(i,:),ep,D);  %ska detta f ens vara med?
    K = assem(edof(i, :), K, Ke);
    [f_b_top, Kce_top] = findRV(e, t, p, a_air, i, 3, T_0);
    [Kc, fb] = assem(edof(i, :), Kc, Kce_top, fb, f_b_top);
end

K = K + Kc;
f = fb;

a = solveq(K, f, bc_lower);        %max(a)=50.1969, hur?!?
ed=extract(edof,a);
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
    Kc = zeros(ndof);
    fb = zeros(ndof,1);

    cur_temp = T_p(tt(j)+dt, T_ref, T_0,t_p);
    bc_lower(:,2) = cur_temp;
    for i = 1:nelm
        
        Ke = flw2te(ex(i,:), ey(i,:), ep, D); 
        K = assem(edof(i, :), K, Ke);
        Ce = plantml(ex(i,:), ey(i,:), cp*rho); %Delar man cp * rho med 10^4 ser det bra ut...    
        C = assem(edof(i,:),C,Ce);
        [f_b_top, Kce_top] = findRV(e, t, p, a_air, i, 3, T_0);
        [Kc, fb] = assem(edof(i, :), Kc, Kce_top, fb, f_b_top);
    end
    K = K + Kc;
    f = fb;
    %g?r divisionen mha solveq ist f?r anew f?r att f? hj?lp med RV
    [a, Q] = solveq(C+K*dt, C*a + f*dt, bc_lower);
    %a = solveq(K*dt, f*dt, bc_lower);      L?ser utan C, blir station?rt
    if(tt(j) == 2.5 ||tt(j) == 5 ||tt(j) == 7.5) 
        figure
        ed=extract(edof,a);
        fill(ex',ey',ed');
        colorbar;
    end
end
