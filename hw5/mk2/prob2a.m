clear;clc;

L = 6;     H = 2;  t = 0.1;

E = 75e6;   nu = 0.3;   P0 = -100e3;

Nx = 60;    Ny = 20;

p = 3;

rho = ones(Nx*Ny,1);

fd = forward_diff_Pavg(rho);

edof = list_dofs(1:Nx*Ny,Nx);

[ d , K , F , ~ , free_eq ] = fem_hw5( L , H  , Nx , Ny , edof , t , rho , p , E , nu , P0 );

adj = adj_method_Pavg(rho);