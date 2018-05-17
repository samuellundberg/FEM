%% geometri & konstanter
mesh = load('mesh1.mat');                        %Laddar den sparade genomerin av en rektangel
e = mesh.e;                                      %Randen        (edges)
p = mesh.p;                                      %Noderna       (points)
t = mesh.t;                                      %Elementen     (triangels)
[m,n] = size(p);

%Allokerar matriser f�r randvilkor
bc_lower = [];
bc_upper = [];
bc_left = [];
bc_right = [];

%Skapar matriser f�r randvilkor som kommer avn�ndas senare
for i = 1:n
    if p(2,i) == -1.5         
        bc_lower = [bc_lower; i, 0];
    end
    if p(2,i) == 1.5
        bc_upper = [bc_upper; i, 0];
    end
    if p(1,i) == -15.0
        bc_left = [bc_left; i, 0];
    end
    if p(1,i) == 15.0
        bc_right = [bc_right; i, 0];
    end
end

%Allokerar matris f�r formfunktioner som kommer anv�ndas f�r att ber�kna
%f_b och Kc
fbvec_top = zeros(3,length(t(1,:)));                                                 
fbvec_bot = zeros(3,length(t(1,:)));

%Skapar matris med formfunktion f�r element i om den ligger p� �vre randen
for i = 1:length(t(1,:))
    fbvec_top = assemForm(t,p,i,fbvec_top,1);
    fbvec_bot = assemForm(t,p,i,fbvec_top,2);
end

%Boj balken
for i = 1:n
    p(2,i)= p(2,i)+(-6.5*sin(pi/15*p(1,i)));
end

p = p/1000;                                 %G�r om till postionen av koordinater till SI-enhet [m]

%Plockar fram n?dv?ndiga parametrar som vi kommer anv?nda 
nelm=length(t(1,:));                        %Antalet element
edof(:,1)=1:nelm;                           %L?gger in elementen i Edof p? rad 1
edof(:,2:4)=t(1:3,:)';                      %L?gger in elementens noder p? rad 2,3,4
coord=p';                                   %D?per om och transponerar punkterna
ndof=length(p);                             %Ber?knar antalet frihetsgrader
[ex,ey]=coordxtr(edof,coord,(1:ndof)',3);   %Tar fram x- och ykoordinater f?r elementen

%Plottar genometrin
eldraw2(ex,ey,[1 2 0]);

% Parametrar som kommer anv?ndas under projektet
ep = 1;                                      %Tjocklek av geometrin
eq = 0;                                      %V?rmef?rs?rjning per volymenhet
%Begynnelsev?rden och randvilkor
T_0 = 25;                                    %Intiella temperaturen
T_inf = 25;                                  %Temperatpur p� ovansidan av balken
T_ref = 50;                                  %Temperaturen efter tiden t_p = uppstartningstiden
t_p = 10;                                    %Uppstartningstiden
a_air = 1000;                                %Konvektionskoefficient
p_air = 10^10;                               %Lufttryck ovansida
p_ref = 5*10^10;                             %Tryck undersida
%matrialkonstanter
alpha = 1.6*10^(-5);                         %Termisk expansionskonstant
k = 16.2;                                    %Konduktivitetskonstant
D = k*eye(2);                                %Konduktivitetsmatrisen
rho = 7990;                                  %Densitet f?r v?rt material [Kgm^?3]
E = 1.93*10^11;                              %youngs modulus
v = 0.25;                                    %poissons ratio
cp = 500;                                    %specifik v?rmekapicitet
%% Stationary solution
K = zeros(ndof);                             %Allocates K-matrix
f = zeros(ndof,1);                           %Allocates f-vector
bc_lower(:,2) = T_ref;                       %Puts the boundaryconditions at the bottom edge to 50C

%loops over every element and calculates the matrices/vectors Ke, Kce, fb and
%assembelse them into K and f
for i = 1:nelm
    Ke = flw2te(ex(i,:),ey(i,:),ep,D); 
    [fb_top, Kce_top] = fbcKec(p, t, fbvec_top, i, a_air, T_inf);
    Ke = Ke + Kce_top;
    [K, f] = assem(edof(i, :), K, Ke, f, fb_top);
end

figure
statT = solveq(K, f, bc_lower);     
ed=extract(edof,statT);
fill(ex',ey',ed');
colorbar;
%% Numdiff
M = 16;                                      %The number of steping
tend = 10;                                   %End time
dt = tend/M;                                 %Size of time steping
tt = linspace(0,tend, M+1);                  %Time spacing
a = T_0*ones(ndof,1);                        %Allocates initial temperatur

%Loops through the time stepping
for j = 1:M
    K = zeros(ndof);                         %Allocates K-matrix
    C = zeros(ndof);                         %Allocates C-matrix
    f = zeros(ndof,1);                       %Allocates f-vector
    
    %Function that calcutaes T_p = (T_ref-T_0)f(t)+T_0
    %where f(t) = (tt/t_p)^3*(6*(tt/t_p)^2-15*tt/t_p+10);
    cur_temp = T_p(tt(j)+dt, T_ref, T_0,t_p);%Sets the current temperature
    bc_lower(:,2) = cur_temp;                %Sets the current temparture as a boundary condition at the lower bondary
    
    %Loops over every element and calculates the matrices/vectors Ce, Ke, Kce, fb and
    %assembelse them into C, K and f
    for i = 1:nelm
        Ke = flw2te(ex(i,:), ey(i,:), ep, D);       
        Ce = plantml(ex(i,:), ey(i,:), cp*rho);    
        C = assem(edof(i,:),C,Ce);
        [fb_top, Kce_top] = fbcKec(p, t, fbvec_top, i, a_air, T_inf); 
        Ke = Ke + Kce_top;
        [K, f] = assem(edof(i, :), K, Ke, f, fb_top);
    end

    %Solves the implisit time stepping sheme
    a = solveq(C+K*dt, C*a + f*dt, bc_lower);
    %Plots the beam after 2.5, 5 and 7.5 seconds
    if(tt(j) == 2.5 ||tt(j) == 5 ||tt(j) == 7.5)
        figure
        ed=extract(edof,a);
        fill(ex',ey',ed');
        colorbar;
    end
end

%% Stress
ptype = 2;                                   %Plain strain
ep = [ptype, 1];                             %ep = [ptype t], ptype: 1=plain stress, 2=plain strain, thickness t = 1
D_hooke = hooke(ptype, E, v);                %Assembels the D matrix
eq = [0, 0];

bc_left(:,2) = 0;                           %Bonndary conditions at the left side of the beam
bc_right(:,2) = 0;                          %Bonndary conditions at the right side of the beam
bc = [bc_left; bc_right];                   %Assembels a matrix bc with the left and right boundary conditions
K = zeros(2*ndof);                          %Allocates K-matrix
fb = zeros(2*ndof,1);                       %Allocates f-vector
f =  zeros(2*ndof,1);                       %Allocates f0-vector

%Expands edof since every node has 2 degrees of freedom 
c1 = [edof(:,2)*2 - 1, edof(:,2)*2];
c2 = [edof(:,3)*2 - 1, edof(:,3)*2];
c3 = [edof(:,4)*2 - 1, edof(:,4)*2];
edof_wide = [edof(:,1), c1, c2, c3];        %Nanes new edof with 2 degrees of freedom at every node edof_wide

%Loops over every element and calculates the matrices/vectors Ke, fe, fb, the
%mean stationary temperature, calculates epsilon0 and assembelse them into K and f
for i = 1:nelm
    Ke = plante(ex(i,:), ey(i,:), ep, D_hooke);
    dT = (statT(t(1,i))+statT(t(2,i))+statT(t(3,i)))/3 - T_0;
    epsilon0 = alpha*dT*[1 1 0]';
    fe = plantf(ex(i,:), ey(i,:), ep, (D_hooke([1 2 4], [1 2 4])*epsilon0)');
    fb_top = findRV2(e, t, p, p_air, i ,3);%SKRIV OM DENNA FUNKTIONEN S� DEN �R P� SAMMA FORMAT SOM �VER
    fb_bot = findRV2(e, t, p, p_ref, i ,1);
    fe = fe + fb_top + fb_bot;
    [K, f] = assem(edof_wide(i, :), K, Ke, f, fe);
end

u = solveq(K, f, bc);                        %Calculates tension in the beam
ed=extract(edof_wide,u);

%VON MISES

Seff_el = zeros(nelm, 1);                    %Allocates Seff elements vector 
Seff_nod = zeros(ndof, 1);                   %Allocates Seff nodes vector 

%Loops over every element and calculates the matrices/vectors Ke, fe, fb, the
%mean stationary temperature, calculates epsilon0 and assembelse them into K and f
%KOLLA UPP VAD SAKERNA HETER!!!!
for i = 1:nelm
    [es, et] = plants(ex(i,:), ey(i,:), ep, D_hooke([1 2 4], [1 2 4]), ed(i,:));
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