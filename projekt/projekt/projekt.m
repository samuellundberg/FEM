%% geometri & konstanter
mesh = load('mesh1.mat');                        %Loads mesh of a rectangle 
e = mesh.e;                                      %Edgees        (edges)
p = mesh.p;                                      %Nodes         (points)
t = mesh.t;                                      %Elements      (triangels)
[m,n] = size(p);

%Allocates boundary condition matrices
bc_lower = [];
bc_upper = [];
bcLR2dof = [];

%Creates matrise for bondary condidtions that will be used later 
for i = 1:n
    if p(2,i) == -1.5         
        bc_lower = [bc_lower; i, 0];
    end
    if p(2,i) == 1.5
        bc_upper = [bc_upper; i, 0];
    end
    if (p(1,i) == 15 || p(1,i) == -15)
        bcLR2dof = [bcLR2dof; 
                   i*2 - 1, 0;
                   i*2    , 0];
    end
end

%Alocates matrix for form functions that will be used to calculate the bondary 
%vector f_b and stiffness matrix Kc
fbvec_top = zeros(3,length(t(1,:)));                                                 
fbvec_bot = zeros(3,length(t(1,:)));

%Creates matrix with form functions for each element on the upper and lower
%edge 
for i = 1:length(t(1,:))
    fbvec_top = assemForm(t,p,i,fbvec_top,3);
    fbvec_bot = assemForm(t,p,i,fbvec_bot,1);
end

%Creates the geometrically simplified plate
for i = 1:n
    p(2,i)= p(2,i)+(-5*sin(pi/15*p(1,i)));
end

p = p/1000;                                 %Sets the geomety to meter [M]

%Creates FEM parameters
nelm=length(t(1,:));                        %Number of elements
edof(:,1)=1:nelm;                           %Puts the elementents i edof at row 1
edof(:,2:4)=t(1:3,:)';                      %Puts the element nodes at row 2,3,4 in edof
coord=p';                                   %Renames and transpoes the nodes vector
ndof=length(p);                             %Calculates the number of degrees of freedom
[ex,ey]=coordxtr(edof,coord,(1:ndof)',3);   %Calculates the x- and ycoordinates for the elements

%Plots the genometry
eldraw2(ex,ey,[1 2 0],edof(:,1));

% Parameters that will be used
ep = 1;                                      %Thickness of the geometry
eq = 0;                                      %The heat supply per unit volume
%Initial and boundary conditions.
T_0 = 25;                                    %Initial temperature, steel & refrigerant 
T_inf = 25;                                  %Surrounding temperature, air
T_ref = 50;                                  %Operating temperature, refrigerant
t_p = 10;                                    %Start up time 
a_air = 1000;                                %Convection coefficient, air
p_air = 10^7;                                %Pressure, air
p_ref = 5*10^7;                              %Pressure, refrigerant
%Material properties of 316 Stainless steel
alpha = 1.6*10^(-5);                         %Thermal expansion coefficient 
k = 16.2;                                    %Thermal conductivity 
D = k*eye(2);                                %Thermal conductivity matrix
rho = 7990;                                  %Density
E = 1.93*10^11;                              %Young’s modulus
v = 0.25;                                    %Poisson’s ratio 
cp = 500;                                    %Specific heat capacity
%% Stationary solution
K = zeros(ndof);                             %Allocates K-matrix
f = zeros(ndof,1);                           %Allocates f-vector
bc_lower(:,2) = T_ref;                       %Puts the boundaryconditions at the bottom edge to 50C

%loops over every element and calculates the matrices/vectors Ke, Kce, fb and
%assembelse them into K and f
for i = 1:nelm
    Ke = flw2te(ex(i,:),ey(i,:),ep,D); 
    [fb_top, Kce_top] = assemFbcKec(p, t, fbvec_top, i, a_air, T_inf);
    Ke = Ke + Kce_top;
    [K, f] = assem(edof(i, :), K, Ke, f, fb_top);
end

%Plots figur
figure(1)
statT = solveq(K, f, bc_lower);     
edstat=extract(edof,statT);
fill(ex',ey',edstat'); colormap('Jet');
title({'Stationary solution temperature distribution, in celcius [°C]', 'of the geometrically simplified plate.'}); xlabel('Legnth of plate. [Meter]'); ylabel('Hight of plate.[Meter]');
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
        [fb_top, Kce_top] = assemFbcKec(p, t, fbvec_top, i, a_air, T_inf); 
        Ke = Ke + Kce_top;
        [K, f] = assem(edof(i, :), K, Ke, f, fb_top);
    end

    %Solves the implisit time stepping sheme
    a = solveq(C+K*dt, C*a + f*dt, bc_lower);
    %Plots the beam after 2.5, 5 and 7.5 seconds
    if(tt(j) == 2.5 ||tt(j) == 5 ||tt(j) == 7.5)
        figure
        ed=extract(edof,a);
        fill(ex',ey',ed'); colormap('Jet');
        title({'Temperature distribution, in celcius [°C], of the geometrically',' simplified plate at time t= X after start up.'}); xlabel('Legnth of plate. [Meter]'); ylabel('Hight of plate.[Meter]');
        colorbar;
    end
end

%% Stress
ptype = 2;                                   %Plain strain
ep = [ptype, 1];                             %ep = [ptype t], ptype: 1=plain stress, 2=plain strain, thickness t = 1
D_hooke = hooke(ptype, E, v);                %Assembels the D matrix
eq = [0, 0];

bc_left(:,2) = 0;                            %Bonndary conditions at the left side of the beam
bc_right(:,2) = 0;                           %Bonndary conditions at the right side of the beam
bc = [bc_left; bc_right];                    %Assembels a matrix bc with the left and right boundary conditions
K = zeros(2*ndof);                           %Allocates K-matrix
fb = zeros(2*ndof,1);                        %Allocates f-vector
f =  zeros(2*ndof,1);                        %Allocates f0-vector
ex_wide = [ex(:,1), ex(:,1), ex(:,2), ex(:,2), ex(:,3), ex(:,3)];
ey_wide = [ey(:,1), ey(:,1), ey(:,2), ey(:,2), ey(:,3), ey(:,3)];

%Expands edof since every node has 2 degrees of freedom 
c1 = [edof(:,2)*2 - 1, edof(:,2)*2];
c2 = [edof(:,3)*2 - 1, edof(:,3)*2];
c3 = [edof(:,4)*2 - 1, edof(:,4)*2];
edof_wide = [edof(:,1), c1, c2, c3];         %Names new edof with 2 degrees of freedom at every node edof_wide

%Loops over every element and calculates the matrices/vectors Ke, fe, fb, the
%mean stationary temperature, calculates epsilon0 and assembelse them into K and f
for i = 1:nelm
    Ke = plante(ex(i,:), ey(i,:), ep, D_hooke);
    dT = sum(edstat(i,:))/3 - T_0;
    epsilon0 = (1+v)*alpha*dT*[1 1 0]';
    fe = plantf(ex(i,:), ey(i,:), ep, (D_hooke([1 2 4], [1 2 4])*epsilon0)');
    fb_top = assemFbcStrain(p, t, fbvec_top, i, p_air, 3);
    fb_bot = assemFbcStrain(p, t, fbvec_bot, i, p_ref, 1);
    fe = fe + fb_top + fb_bot;
    [K, f] = assem(edof_wide(i, :), K, Ke, f, fe);
end

u = solveq(K, f, bcLR2dof);                  %Calculates tension in the beam
ed=extract(edof_wide,u);                     %Node displacement

%VON MISES
Seff_el = zeros(nelm, 1);                    %Allocates von Mises stress vector for elements 
Seff_nod = zeros(ndof, 1);                   %Allocates von Mises stress vector for nodes

%Loops over every element and calculates the von Mises stress for each
%element
for i = 1:nelm
    [es, et] = plants(ex(i,:), ey(i,:), ep, D_hooke([1 2 4], [1 2 4]), ed(i,:)); %Stresses and strains
    dT = sum(edstat(i,:))/3 - T_0;
    sigma_zz= v*(es(1) + es(2)) - alpha*E*dT; %side 256 i course litterature
    Mises = sqrt(es(1,1)^2 + es(1,2)^2 + sigma_zz^2 - es(1,1)*es(1,2) - es(1,1)*sigma_zz - es(1,2)*sigma_zz + 3*es(1,3)^2);
    Seff_el(i,1) = Mises;
    if i == 43 
        Mises
    end
end

%Sets the von Mises stress for each element of freedom
for i=1:ndof
    [c0,c1]=find(edof(:,2:4)==i);
    Seff_nod(i,1)=sum(Seff_el(c0))/size(c0,1);
end

%Plots the effective von Mises stress field
figure(1)
ed1 = extract(edof, Seff_nod);
fill(ex',ey',ed1')
title({'The effective von Mises stress field, in pascal [Pa] of the geometrically','simplified plate at stationary conditions'}); xlabel('Legnth of plate. [Meter]'); ylabel('Hight of plate.[Meter]');
colorbar; 

%Plots the displacement field
figure(2)
eldraw2(ex, ey, [1 2 1]);
eldisp2(ex,ey, ed, [1 4 1], 10);
title({'The displacement field of the of the geometrically simplified plate', 'amplified a factor 10 at stationary conditions'}); xlabel('Legnth of plate. [Meter]'); ylabel('Hight of plate.[Meter]');