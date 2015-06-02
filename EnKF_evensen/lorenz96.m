% Lorenz-96 model

% dx(i)/dt = x(i-1)*(x(i+1)-x(i-2))-x(i)+F
% x(-1)=x(Nx-1)
% x(0)=x(Nx)
% x(1)=x(Nx+1)

% Nx=40;
% F=8.0

function yprime = lorenz96(~, y)

N = 40;  %  System states
F = 8.0 * ones(N, 1);  % Chaotic system parameter

% F=4.0*ones(N,1);
A = (-1) * eye(N, N);

% cyclic boundary condition
A(N, 1) = y(N - 1);
A(1, N - 1) = -y(N);
A(2, N) = - y(1);
A(1, 2) = y(N);


for i = 2 : N - 1
   A(i, i + 1) = y(i - 1);
end

for i = 3 : N
   A(i, i - 2) = - y(i - 1); 
end

% Structure of A
% -1      y(N)          ...   -y(N)      0
%  0       -1    y(2-1)  ...            -y(1)
% -y(i-1) ...    -1    y(i-1) ... 
% .
% .
% .
% .            -y(i-1) -1 ...          y(i-1)
% y(N-1)  ...               -y(i-1)  0   -1



yprime = A * y + F;



