function dfdrhoi = adj_method_Pavg(r)
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

    edof = list_dofs(1:Nx*Ny,Nx);
    
    %       Solve the FEM problem for the initial design point
    [ d , K , ~ , fix_eq , free_eq ] = fem_hw5( xmax , ymax  , Nx , Ny , edof , t , r , p , E , nu , P0 );

    %       Adjoint vector
    sy = zeros(dofs,1);
    
    %       Explicit dep. of func of interest on state vars (disp.) dpfdpu
    L = zeros(dofs,1);
    Fx_n = 2*(1 + (Nx+1)*( (1:Ny+1) - 1))-1;    N = length(Fx_n);
    L(Fx_n) = 1;
    
    %       Prescribed part of adjoint vector
    sy(fix_eq) = L(fix_eq)/N;
    
    %       Free part of adjoint vector
    sy(free_eq) = -transpose(K(free_eq,free_eq))\transpose(K(fix_eq,free_eq))*sy(fix_eq);
    
    dfdrhoi = zeros(Nx*Ny,1);
    elemK = elementK( xmax , ymax , Nx , Ny , t , E , nu);
    for e = 1:Nx*Ny
        eqn_num = edof(e,:);
        dfdrhoi(e) = transpose(sy(eqn_num))*(p*r(e)^(p-1))*elemK*d(eqn_num);
    end

end

% %--------------------------------------------------------------------------
% %   Material constants (in SI units)
% %--------------------------------------------------------------------------
% E = 75e6;   nu = 0.3;   t = 0.1;    p = 3;
% P0 = -100e3;
% %--------------------------------------------------------------------------
% %   Mesh params
% %--------------------------------------------------------------------------
% xmax = 6;  ymax = 2;
% Nx = 60;   Ny = 20;
% %--------------------------------------------------------------------------
% 
% dofs = 2 * ( Nx + 1 ) * ( Ny + 1 );
% r = 1*ones(Nx*Ny,1);
% edof = list_dofs(1:Nx*Ny,Nx);
% tic
% 
% %       Solve the FEM problem for the initial design point
% [ d , K , F , fix_eq , free_eq ] = fem_hw5( xmax , ymax  , Nx , Ny , edof , t , r , p , E , nu , P0 );
% toc
% tic
% %       Adjoint vector
% sy = zeros(dofs,1);
% 
% %       Explicit dep. of func of interest on state vars (disp.) dpfdpu
% L = zeros(dofs,1);
% Fx_n = 2*(1 + (Nx+1)*( (1:Ny+1) - 1))-1;    N = length(Fx_n);
% L(Fx_n) = 1;
% 
% %       Prescribed part of adjoint vector
% sy(fix_eq) = L(fix_eq)/N;
% 
% %       Free part of adjoint vector
% sy(free_eq) = -transpose(K(free_eq,free_eq))\transpose(K(fix_eq,free_eq))*sy(fix_eq);
% 
% dfdrhoi = zeros(Nx*Ny,1);
% elemK = elementK( xmax , ymax , Nx , Ny , t , E , nu);
% for e = 1:Nx*Ny
%     eqn_num = edof(e,:);
%     dfdrhoi(e) = transpose(sy(eqn_num))*(p*r(e)^(p-1))*elemK*d(eqn_num);
% end
% 
% toc
