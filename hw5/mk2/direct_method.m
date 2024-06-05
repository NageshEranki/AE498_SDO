clear;clc;

%   Define constants
%--------------------------------------------------------------------------
Nx = 5;     Ny = 2;

L = 0.5;    H = 0.2;    t = 0.05;

E = 150e9;  nu = 0.3;   P0 = -2e9;
%--------------------------------------------------------------------------

dofs = 2*(Nx+1)*(Ny+1);

tic

%   Solve the FEM prob at the current design point
[ d , K , F , fix_eq , free_eq ] = fem_pstress( L , H , Nx , Ny , t , E , nu , P0 );

%   Explicit dependence of residual (R) on the state vars (displacements)
dRdy = zeros(dofs,dofs);
dRdy(fix_eq,fix_eq) = -eye(length(fix_eq),length(fix_eq));
dRdy(free_eq,free_eq) = K(free_eq,free_eq);
dRdy(fix_eq,free_eq) = K(fix_eq,free_eq);


%   Explicit dependence of residual (R) on the design vars (thicknesses)
dRdti = zeros(dofs,Nx*Ny);
for e = 1:Nx*Ny
    dKdti = zeros(dofs,dofs);
    eqn_num = list_dofs( e , Nx );
    elemK = elementK( L , H , Nx , Ny , 1 , E , nu );
    dKdti(eqn_num,eqn_num) = dKdti(eqn_num,eqn_num) + elemK;
%     ti = t*zeros(Nx*Ny,1);
%     ti(e) = 1;
%     dKdti = globalK( L , H , Nx , Ny , ti , E , nu );
    dRdti(:,e) = dKdti*d;
end

%   Implicit dependence of state vars (displacements) on design vars
%   (thicknesses)
dydti = -dRdy\dRdti;

%   Explicit dependence of func of interest on state vars (displacements)
dfdy = zeros(1,dofs);
dfdy(2*Nx+2) = 1;

%   Final sensitivity formula
dfdti = dfdy*dydti;

toc