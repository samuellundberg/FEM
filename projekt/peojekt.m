%% geometrin
mesh = load('mesh1.mat');    %Laddar den sparade genomerin av en rektangel
e = mesh.e;                  %Randen        (edges)
p = mesh.p;                  %Noderna       (points)
t = mesh.t;                  %Elementen     (triangels)
[m,n] = size(p);

%Flyttar alla noder med avståndet -6.5*sin(pi*x/15) för att få geometrin av
%den önskade sinusvågen sinusvåg
for i = 1:n
    p(2,i)= p(2,i)+(-6.5*sin(pi/15*p(1,i)));
end 

%Plockar fram nödvändiga parametrar som vi kommer använda 

nelm=length(t(1,:));                        %Antalet element
edof(:,1)=1:nelm;                           %Lägger in elementen i Edof på rad 1
edof(:,2:4)=t(1:3,:)';                      %Lägger in elementens noder på rad 2,3,4
coord=p' ;                                  %Döper om och transponerar punkterna
ndof=max(max(t(1:3,:)));                    %Beräknar antalet frihetsgrader
[ex,ey]=coordxtr(edof,coord,(1:ndof)',3);   %Tar fram x- och ykoordinater för elementen

%Plottar genometrin
eldraw2(ex,ey,[1,2,0], edof(:,1)) ;                    
%% Parametrar som kommer användas under projektet
ep = 1;                 %Tjocklek av geometrin
eq = 1;                 %Värmeförsörjning per volymenhet

T_0=25;                 %Intiella temperaturen
T_inf = 25;             %Temperaturen "långt bort"
T_ref = 50;             %Temperaturen efter tiden t_p = uppstartningstiden

t_p = 10;               %Uppstartningstiden
    
%Funktion som beskriver temperaturen längs rand efter tiden t
%t = 1
%f = @(t) (t/t_p).^3.*(6.*(t/t_p).^2-15.*t/t_p+10);
%T_p = (T_ref-T_0)*f(t)+T_0;


a_air = 1000;           %Konvektionskoefficient 
%q_n = a_air*(t-T_inf); %Newtons konvektionslag, beskriver värmeflöde
alpha = 1.6e-5;         %Termisk expansionskonstant
k = 16.2;               %Konduktivitetskonstant
D = k*eye(2);           %Konduktivitetsmatrisen
rho = 7990;             %Densitet för vårt material [Kgm^?3]



%%

K = zeros(ndof);
f = zeros(ndof,1);
Kce = zeros(3);
for i = 1:nelm
   [Ke, fe] = flw2te(ex(i,:),ey(i,:),ep,D,eq);
   Kce = assemKce(ex(i,:),ey(i,:),alpha);
   
   Ce=plantml(ex(i,:),ey(i,:),rho);                 %Skapar Ce-matrisen från kurshemsidan
    
    [K, f] = assem(edof(i, :), K, Ke, f, fe);
end

a = solveq(K,f);

ed=extract(edof,a);
colormap(hot);
fill(ex',ey',ed');
colorbar;
