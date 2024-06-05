function dfdrhoi = adj_method_Ucy(rho,d,K,free_eq,xmax,ymax,Nx,Ny,t,E,nu,p)
    
    dofs = 2 * ( Nx + 1 ) * ( Ny + 1 );

    %       Adjoint vector
    sy = zeros(dofs,1);

    
    %       Explicit dep. of func of interest on state vars (disp.) dpfdpu

    L = sparse(dofs,1);
    L(end) = 1;
    
    sy(free_eq) = -transpose(K(free_eq,free_eq))\L(free_eq);
    elemK = elementK( xmax , ymax , Nx , Ny , t , E , nu );
    dfdrhoi = zeros( Nx*Ny , 1 );
    for e = 1:Nx*Ny    
    
        eqn_num = list_dofs( e , Nx );
    
        dfdrhoi(e) = transpose(sy(eqn_num)) * (p*rho(e)^(p-1) * elemK )  * d(eqn_num);
    
    end
end

%{
clear;clc;
%--------------------------------------------------------------------------
%   Material constants (in SI units)
%--------------------------------------------------------------------------
E = 75e6;   nu = 0.3;   t = 0.1;    p = 3;
P0 = -100e3;
%--------------------------------------------------------------------------
%   Mesh params
%--------------------------------------------------------------------------
xmax = 6;  ymax = 2;
Nx = 60;   Ny = 20;
%--------------------------------------------------------------------------

dofs = 2 * ( Nx + 1 ) * ( Ny + 1 );
r = 0.4*ones(Nx*Ny,1);

tic

%       Solve the FEM problem for the initial design point
edof=list_dofs(1:Nx*Ny,Nx);
[ d , K , ~ , ~ , free_eq ] = fem_hw5( xmax , ymax , Nx , Ny , edof , t , r , p , E , nu , P0 );
toc
tic
%       Adjoint vector
sy = zeros(dofs,1);

%       Explicit dep. of func of interest on state vars (disp.) dpfdpu
L = zeros(dofs,1);
L(end) = 1;

%       Extract those entries corresponding to free dofs
%dfduff = L(free_eq);

%       Extract from K the entries corresp. to free dofs
%Kff = K(free_eq,free_eq);

sy(free_eq) = -transpose(K(free_eq,free_eq))\L(free_eq);
toc
tic
elemK = elementK( xmax , ymax , Nx , Ny , t , E , nu );
dfdrhoi = zeros( Nx*Ny , 1 );
for e = 1:Nx*Ny    

    eqn_num = list_dofs( e , Nx );

    dfdrhoi(e) = transpose(sy(eqn_num)) * (p*r(e)^(p-1) * elemK )  * d(eqn_num);

end

toc
%}
