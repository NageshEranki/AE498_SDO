function fd = forward_diff(h)
% clear;clc;

%   Define constants
%--------------------------------------------------------------------------
Nx = 5;     Ny = 2;

% h = 1e-5*1e-2;

L = 0.5;    H = 0.2;    t = 0.05*ones(Nx*Ny) + h * eye(Nx*Ny);

E = 150e9;  nu = 0.3;   P0 = -2e9;


%--------------------------------------------------------------------------
tic

fd = zeros(Nx*Ny,1);
d = fem_pstress( L , H , Nx , Ny , 0.05 , E , nu , P0 );
d0 = d(2*Nx+2);
for i = 1:Nx*Ny
    d = fem_pstress( L , H , Nx , Ny , t(:,i) , E , nu , P0 );
    di = d(2*Nx+2);
    fd(i) = (di-d0)/h;
end

% fd;

% toc

end