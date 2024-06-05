clear;clc;

%   Define constants
%--------------------------------------------------------------------------
Nx = 5;     Ny = 2;

L = 0.5;    H = 0.2;    t = 0.05;

E = 150e9;  nu = 0.3;   P0 = -2e9;

%--------------------------------------------------------------------------

dofs = 2 * ( Nx + 1 ) * ( Ny + 1 );

tic

%       Solve the FEM problem for the initial design point
[ d , K , ~ , ~ , free_eq ] = fem_pstress( L , H , Nx , Ny , t , E , nu , P0 );

%       Adjoint vector
sy = zeros(dofs,1);

%       Explicit dep. of func of interest on state vars (disp.)
dfdu = zeros(dofs,1);
dfdu(2*Nx+2) = 1;

%       Extract those entries corresponding to free dofs
dfduff = dfdu(free_eq);

%       Extract from K the entries corresp. to free dofs
Kff = K(free_eq,free_eq);

sy(free_eq) = -transpose(Kff)\dfduff;

dfdti = zeros( 1 , Nx*Ny );
for e = 1:Nx*Ny

    dKdti = zeros(dofs,dofs);

    elemK = elementK( L , H , Nx , Ny , 1 , E , nu );

    eqn_num = list_dofs( e , Nx );

    dKdti(eqn_num,eqn_num) = elemK;

    dfdti(e) = transpose(sy) * dKdti * d;

end

toc