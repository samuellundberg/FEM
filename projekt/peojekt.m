%% geometrin
mesh = load('mesh1.mat');    %Laddar den sparade genomerin av en rektangel
e = mesh.e;                  %Randen        (edges)
p = mesh.p;                  %Noderna       (points)
t = mesh.t;                  %Elementen     (triangels)
[m,n] = size(p);

%Flyttar alla noder med avst�ndet -6.5*sin(pi*x/15) f�r att f� geometrin av
%den �nskade sinusv�gen sinusv�g
for i = 1:n
    p(2,i)= p(2,i)+(-6.5*sin(pi/15*p(1,i)));
end 

%Plockar fram n�dv�ndiga parametrar som vi kommer anv�nda 

nelm=length(t(1,:));                        %Antalet element
edof(:,1)=1:nelm;                           %L�gger in elementen i Edof p� rad 1
edof(:,2:4)=t(1:3,:)';                      %L�gger in elementens noder p� rad 2,3,4
coord=p' ;                                  %D�per om och transponerar punkterna
ndof=max(max(t(1:3,:)));                    %Ber�knar antalet frihetsgrader
[ex,ey]=coordxtr(edof,coord,(1:ndof)',3);   %Tar fram x- och ykoordinater f�r elementen

%Plottar genometrin
eldraw2(ex,ey,[1,2,0], edof(:,1)) ;                    
%% Parametrar som kommer anv�ndas under projektet
ep = 1;                 %Tjocklek av geometrin
eq = 1;                 %V�rmef�rs�rjning per volymenhet

T_0=25;                 %Intiella temperaturen
T_inf = 25;             %Temperaturen "l�ngt bort"
T_ref = 50;             %Temperaturen efter tiden t_p = uppstartningstiden

t_p = 10;               %Uppstartningstiden
    
%Funktion som beskriver temperaturen l�ngs rand efter tiden t
%t = 1
%f = @(t) (t/t_p).^3.*(6.*(t/t_p).^2-15.*t/t_p+10);
%T_p = (T_ref-T_0)*f(t)+T_0;


a_air = 1000;           %Konvektionskoefficient 
%q_n = a_air*(t-T_inf); %Newtons konvektionslag, beskriver v�rmefl�de
alpha = 1.6e-5;         %Termisk expansionskonstant
k = 16.2;               %Konduktivitetskonstant
D = k*eye(2);           %Konduktivitetsmatrisen
rho = 7990;             %Densitet f�r v�rt material [Kgm^?3]



%%

K = zeros(ndof);
f = zeros(ndof,1);
Kce = zeros(3);
for i = 1:nelm
   [Ke, fe] = flw2te(ex(i,:),ey(i,:),ep,D,eq);
   Kce = assemKce(ex(i,:),ey(i,:),alpha);
   Ke = Ke + Kce;
   Ce=plantml(ex(i,:),ey(i,:),rho);                 %Skapar Ce-matrisen fr�n kurshemsidan
    
   [K, f] = assem(edof(i, :), K, Ke, f, fe);
end

a = solveq(K,f);

ed=extract(edof,a);
colormap(hot);
fill(ex',ey',ed');
colorbar;
%% test
ex = [0 2 0];
ey = [0 1 1];
alpha = 10;
res = assemKce(ex',ey',alpha)