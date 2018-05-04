%%
n = 20;
x = linspace(2,8,n+1);
L = (x(end)-x(1))/n;
A = 10; k = 5; q=15; Q=100;
Ke = A*k/L*[1 -1 ; -1 1];
edof=zeros(n,3);
for i=1:n
    edof(i,:) = [i i i+1];
end
K = zeros(n+1);
f = zeros(n+1,1); f(end) = -q*A;
fe = [Q*L/2 ; Q*L/2];
for i = 1:n
    [K, f] = assem(edof(i, :), K, Ke, f, fe);
end
bc = [1 0];
T = solveq(K,f,bc);

%%
ep = 1;
eq = 1;
k = 1;
D = k*eye(2);
K = zeros(ndof);
f = zeros(ndof,1);

for i = 1:nelm
   [Ke, fe] = flw2te(ex(i,:),ey(i,:),ep,D,eq);
   [K, f] = assem(edof(i, :), K, Ke, f, fe);
end

bc=[bc3; bc4; bcarc];
T = solveq(K,f,bc);

ed=extract(edof,T);
colormap(hot);
fill(ex',ey',ed');
colorbar;