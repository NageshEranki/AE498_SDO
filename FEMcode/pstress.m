clear;clc;


%   Define constants
%--------------------------------------------------------------------------
Nx = 25;     Ny = 10;

L = 0.5;    H = 0.2;    t = 0.05*ones(Nx*Ny);

E = 150e9;  nu = 0.3;   P0 = -2e9;
%--------------------------------------------------------------------------


%   Total number of degrees of freedom
dofs = 2*(Nx+1)*(Ny+1);


%       Assemble global K
%--------------------------------------------------------------------------
K = globalK( L , H , Nx , Ny , t , E , nu );
%--------------------------------------------------------------------------


%   Global force vector
%--------------------------------------------------------------------------
F = globalF( L , H , Nx , Ny , 0.05 , P0 );
%--------------------------------------------------------------------------


%   List nodes which are fixed
%--------------------------------------------------------------------------
fix_n = 1 + (Nx+1)*( (1:Ny+1) - 1);
fix_eq = zeros(1,2*length(fix_n));
for i = 1:length(fix_n)
    fix_eq(2*i-1) = 2*fix_n(i)-1; fix_eq(2*i) = 2*fix_n(i);
end
%--------------------------------------------------------------------------


%   List eqn numbers corresp. to free dofs
free_eq=setdiff(1:dofs,fix_eq);

%   extract rows and columns of the free dofs
Kred=K(free_eq,free_eq);
Fred=F(free_eq);

dred=Kred\Fred;


%   Rebuild d
d=zeros(dofs,1);
for i=1:length(free_eq)
    d(free_eq(i))=dred(i);
end

straine = 0.5*d'*K*d

%   Plot results
[xi,yi] = ndgrid(linspace(0,L,Nx+1),linspace(0,H,Ny+1));
xi = xi(:); yi = yi(:);
% d2=d;
d2 = (0.1*L)/max(abs(d),[],'all')*d;

dx = d2(1:2:end);    dy = d2(2:2:end);
plot(xi+dx,yi+dy,'xb')
axis([-0.2 L+0.2 -0.2 0.2+H])