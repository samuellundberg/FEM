function Kce=assemKce(ex,ey,a_air)

%-------------------------------------------------------------
% 
%  Beräknar: Kce=a_air*int(N^T*N)dS
%
% INPUT:  ex,ey;       Element coordinates
%	
%	  a_air
%
% OUTPUT: Kce :      Matix 3 x 3
%-------------------------------------------------------------

%Beräknar längden för ett element 
Legnth=sqrt((ex(1)-ex(2))^2+(ey(1)-ey(2))^2) + ...
       sqrt((ex(1)-ex(3))^2+(ey(1)-ey(3))^2) + ...
       sqrt((ex(2)-ex(3))^2+(ey(2)-ey(3))^2);


L1=[0.5 0 0.5];
L2=[0.5 0.5 0];
L3=[0 0.5 0.5];

NtN=zeros(6);


for i=1:3
	NtN=NtN+1/3*[L1(i)^2 0 L1(i)*L2(i) 0 L1(i)*L2(i) 0
		 		0 L1(i)^2 0 L1(i)*L2(i) 0 L1(i)*L2(i)
		  		L2(i)*L1(i) 0 L2(i)^2 0 L2(i)*L3(i) 0
		  		0 L2(i)*L1(i) 0 L2(i)^2 0 L2(i)*L3(i)
		 		L3(i)*L1(i) 0 L3(i)*L2(i) 0 L3(i)^2 0
		  		0 L3(i)*L1(i) 0 L3(i)*L2(i) 0 L3(i)^2];
end

Me1=NtN*Legnth*a_air;


Kce=Me1([1 3 5],[1 3 5]);       

