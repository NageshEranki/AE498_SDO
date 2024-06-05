function fd = forward_diff_Pavg(r0)
    %--------------------------------------------------------------------------
    %   Material constants (in SI units)
    %--------------------------------------------------------------------------
    E = 75e6;   nu = 0.3;   t = 0.1;    p = 3;
    P0 = -100e3;
    %--------------------------------------------------------------------------
    %   Mesh params
    %--------------------------------------------------------------------------
    L = 6;  H = 2;
    Nx = 60;   Ny = 20;
    %--------------------------------------------------------------------------
    %   Elements densities
    edof=list_dofs(1:Nx*Ny,Nx);
    %   Step-size
    h = 1e-7;
    %--------------------------------------------------------------------------
    
    fd = zeros(Nx*Ny,1);
    [ d , K , ~ , fix_eq , free_eq ] = fem_hw5( L , H  , Nx , Ny , edof , t , r0 , p , E , nu , P0 );
    %   DOUBLE CHECK! Dofs corresponding to Fx on left surface only
    Fp = zeros(2*(Nx+1)*(Ny+1));
    Fp(fix_eq) = K(fix_eq,free_eq)*d(free_eq);
    Fx_dof = 2*(1 + (Nx+1)*( (1:Ny+1) - 1))-1;
    f0 = mean(Fp(Fx_dof));
    for i=1:Nx*Ny
        r = r0;
        r(i) = r(i)+h;
        [ di , Ki , ~ , fix_eqi , free_eqi ] = fem_hw5( L , H  , Nx , Ny , edof , t , r , p , E , nu , P0 );
        %   DOUBLE CHECK!
        Fpi = zeros(2*(Nx+1)*(Ny+1));
        Fpi(fix_eq) = Ki(fix_eqi,free_eqi)*di(free_eqi);
        fi = mean(Fpi(Fx_dof));
        fd(i) = (fi-f0)/h;
    end

end
% clear;clc;
% %--------------------------------------------------------------------------
% %   Material constants (in SI units)
% %--------------------------------------------------------------------------
% E = 75e6;   nu = 0.3;   t = 0.1; p = 3;
% P0 = -100e3;
% %--------------------------------------------------------------------------
% %   Mesh params
% %--------------------------------------------------------------------------
% L = 6;  H = 2;
% Nx = 60;   Ny = 20;
% %--------------------------------------------------------------------------
% %   Elements densities
% r0 = 0.4*ones(Nx*Ny,1);
% edof=list_dofs(1:Nx*Ny,Nx);
% %   Step-size
% h = 1e-7*1e-2;
% %--------------------------------------------------------------------------
% tic
% 
% fd = zeros(Nx*Ny,1);
% [ d , K , F , fix_eq , free_eq ] = fem_hw5( L , H  , Nx , Ny , edof , t , r0 , p , E , nu , P0 );
% toc
% %   DOUBLE CHECK! Dofs corresponding to Fx on left surface only
% Fp = zeros(2*(Nx+1)*(Ny+1));
% Fp(fix_eq) = K(fix_eq,free_eq)*d(free_eq);
% Fx_dof = 2*(1 + (Nx+1)*( (1:Ny+1) - 1))-1;
% f0 = mean(Fp(Fx_dof));
% %for i = 1:Nx*Ny
% for i=1:20
% 	tic
%     r = r0;
%     r(i) = r(i)+h;
%     [ di , Ki , Fi , fix_eqi , free_eqi ] = fem_hw5( L , H  , Nx , Ny , edof , t , r , p , E , nu , P0 );
%     %   DOUBLE CHECK!
%     Fpi = zeros(2*(Nx+1)*(Ny+1));
%     Fpi(fix_eq) = Ki(fix_eqi,free_eqi)*di(free_eqi);
%     fi = mean(Fpi(Fx_dof));
%     fd(i) = (fi-f0)/h;
% 	toc
% end
% 
% toc
