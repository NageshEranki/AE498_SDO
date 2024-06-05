clear;clc
%--------------------------------------------------------------------------
%   Material props
%--------------------------------------------------------------------------
E = 75e6;   nu = 0.3;   t = 0.1;    p = 3;
P0 = 100e3;
%--------------------------------------------------------------------------
%   Mesh params
%--------------------------------------------------------------------------
L = 6;  H = 2;
Nx = 60;   Ny = 20;
%--------------------------------------------------------------------------
dofs = 2 * ( Nx + 1 ) * ( Ny + 1 );
edof=list_dofs(1:Nx*Ny,Nx);
%--------------------------------------------------------------------------
%   Intial density
%--------------------------------------------------------------------------
rho = 0.1*ones(Nx*Ny,1);
%--------------------------------------------------------------------------
%   Density filter radius and volume fraction
%--------------------------------------------------------------------------
ep = 0.55;   
vf = 0.4;
%   Filter matrix; w
w = topfilter(L,H,Nx,Ny,ep);
% [x,y] = ndgrid( linspace(0,L,Nx+1) , linspace(0,H,Ny+1));
% w = density_filter(Nx, Ny, 1.5 , x(:) , y(:) , list_nodes(1:Nx*Ny,Nx));
%--------------------------------------------------------------------------
%   Objective and its grad
%--------------------------------------------------------------------------
f = @(rho) objectiveFcn(rho,w,L,H,Nx,Ny,edof,t,p,E,nu,P0);
constr = @(rho) constrFcn(rho,w,L,H,t,Nx,Ny,vf);
%--------------------------------------------------------------------------
Ve = (L*H*t)/(Nx*Ny);   %   Element volume
% A = Ve*ones(1,Nx*Ny);   %   Sum of individual volumes
% b = 0.4*(L*H*t);        %   ''  '' 
A = [];
b = [];
Aeq=[];
beq=[];
lb=1e-5*zeros(1,Nx*Ny);  %   Lower bound on element densities
ub=ones(1,Nx*Ny);       %   Upper bound on element densities
options=optimoptions("fmincon", ...
    "Algorithm","interior-point", ... 
    "SpecifyObjectiveGradient",true, ... 
    "CheckGradients",false, ...
    "FiniteDifferenceType",'central', ...
    "SpecifyConstraintGradient",true, ... 
    "Display","iter", ...
    "StepTolerance",1e-12, ... 
    "MaxIterations",2000, ...
    "MaxFunctionEvaluations",Inf, ...
    "OptimalityTolerance",1e-6, ...
    "ScaleProblem", "obj-and-constr", ...
    "PlotFcn",["optimplotfval","optimplotfirstorderopt"]);
    
x = fmincon(f,rho,A,b,Aeq,beq,lb,ub,constr,options);
show_dens(w*x,L,H,Nx,Ny)

% title("Nx="+Nx+" Ny="+Ny+" eps="+ep+" p="+p+" vf="+vf)
% saveas(gcf,'vol_cont/5.png')

function [f,g] = objectiveFcn(rho,w,L,H,Nx,Ny,edof,t,p,E,nu,P0)
%     show_struct(Nx,Ny,rho);
    %--------------------------------------------------------------------------
    %   Apply density filter
    %--------------------------------------------------------------------------
    rho=w*rho;
    %--------------------------------------------------------------------------
    %   Solve FEM problem based on filtered densities
    %--------------------------------------------------------------------------
    [ d , K , ~ , ~ , free_eq ] = fem_hw5( L , H  , Nx , Ny , edof , t , rho , p , E , nu , P0 );
    %--------------------------------------------------------------------------
    %   Extract Ucy
    %--------------------------------------------------------------------------
    f = d(end);
    %--------------------------------------------------------------------------
    if nargout > 1 
        %------------------------------------------------------------------
        g = adj_method_Ucy(rho,d,K,free_eq,L,H,Nx,Ny,t,E,nu,p);
        %--------------------------------------------------------------------------
        %   Reverse filtering op
        %--------------------------------------------------------------------------
         g = w*g;
    end
end

function [c,ceq,gc,gceq] = constrFcn(rho,w,L,H,t,Nx,Ny,vf)
    V = L*H*t;      %   Total volume of the structure
    Ve = V/(Nx*Ny); %   Volume of each element
    
    %----------------------------------------------------------------------
    %   Filter before evaluating constraints
    %----------------------------------------------------------------------
    rho = w*rho;
    c = Ve*sum(rho) - vf*V;        %   Evaluate constraint
    %----------------------------------------------------------------------
    if nargout > 1
        gc = Ve*ones(Nx*Ny,1);  %   Constraint grad
        gc = w*gc;              %   Reverse filtering op
    end
    ceq=[];
    gceq=[];
end
