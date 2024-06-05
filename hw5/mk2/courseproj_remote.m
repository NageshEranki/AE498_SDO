clear;clc;
tic
str=problem;
esns=connect(str);
%   Define element densities
rho = 0.4*ones(str.Nx*str.Ny,1);
rho(esns.Eg) = 1e-3;
%rho=readmatrix('test.dat');

ep = 1.5*min([str.hx, str.hy]);
p = 3;

w = topfilter(str.L,str.H,str.Nx,str.Ny,ep);
%--------------------------------------------------------------------------
%   Objective and its grad
%--------------------------------------------------------------------------
f = @(rho) objFcn(rho,w,p,str,esns);
constr = @(rho) constrFcn(rho,w,p,str,esns);

%--------------------------------------------------------------------------
A = []; b = []; Aeq=[]; beq=[];
lb=1e-3*ones(1,str.Nx*str.Ny);  %   Lower bound on element densities
ub=ones(1,str.Nx*str.Ny);       %   Upper bound on element densities
ub(esns.Eg) = 1.01e-3;
options=optimoptions("fmincon", ...
    "Algorithm","interior-point", ... 
    "SpecifyObjectiveGradient",true, ... 
    "CheckGradients",false, ...
    "FiniteDifferenceType","central", ...
    "SpecifyConstraintGradient",true, ... 
    "Display","iter", ...
    "MaxIterations",2000, ...
    "MaxFunctionEvaluations",Inf, ...
    "OptimalityTolerance",1e-6, ...
    "StepTolerance",1e-16, ...
    "ScaleProblem","obj-and-constr", ...
    "PlotFcn",{@optimplotfval,@optimplotfirstorderopt});
    
x = fmincon(f,rho,A,b,Aeq,beq,lb,ub,constr,options);
show_struct(str.Nx,str.Ny,w*x)
con=@(rho)constrFcn(rho,w,p,str,esns);
fem=@(rho)femsolve(rho,p,str,esns);
% [d,K,F]=femsolve(rho,p,str,esns);
% temp = sparse(str.dofs,1);
% %   Compute reaction forces
% temp(esns.fix_eq) = K(esns.fix_eq,esns.free_eq)*d(esns.free_eq);
%   Gripping force is the negative
%   of the reaction at the tip
%   of the gripper surface
% f = -temp(esns.dofsg(end))
%   Plot results
% [xi,yi] = ndgrid(linspace(0,str.L,str.Nx+1),linspace(0,str.H,str.Ny+1));
% xi = xi(:); yi = yi(:);
% d2 = (0.1*(str.L))/max(abs(d),[],'all')*d;
% d2 = d;
% 
% dx = d2(1:2:end);    dy = d2(2:2:end);
% plot(xi+dx,yi+dy,'xb')
% axis equal
% [d,K,F] = femsolve(rho,p,str,esns);
% g=F(esns.free_eq)'*d(esns.free_eq)

function [d,K,F] = femsolve(rho,p,str,esns)
    K = globalK( str.L , str.H , str.Nx , str.Ny , ...
        esns.edof , str.T , ...
        rho , p , ...
        str.E , str.nu );
    F = globalF_proj( str.L , str.H , str.Nx , str.Ny , ...
         str.T , str.P0 );
%     F = sparse(str.dofs,1);
%     F(2*(1+str.Ny*(str.Nx+1))) = -1000;
    d=zeros(str.dofs,1);
    d(esns.free_eq) = K(esns.free_eq,esns.free_eq)\F(esns.free_eq);
end

function [f,g] = objFcn(rho,w,p,str,esns)
%      tic
    %--------------------------------------------------------------------------
    %   Apply density filter
    %--------------------------------------------------------------------------
    rho = w*rho;
    %--------------------------------------------------------------------------
    %   Solve FEM problem w filtered densities
    %--------------------------------------------------------------------------
    [d,K,F]=femsolve(rho,p,str,esns);
    %--------------------------------------------------------------------------
    %   Extract force at tip of gripping surface
    %--------------------------------------------------------------------------
    temp = sparse(str.dofs,1);
    %   Compute reaction forces
    temp(esns.fix_eq) = K(esns.fix_eq,esns.free_eq)*d(esns.free_eq);
    %   Gripping force is the negative
    %   of the reaction at the tip
    %   of the gripper surface
    f = -temp(esns.dofsg(end))/abs(sum(F));
%     toc
    %--------------------------------------------------------------------------
    %   Adjoint method for sensitivity
    %--------------------------------------------------------------------------
    if nargout > 1
        g = zeros(str.Nx*str.Ny,1);
        L = sparse(str.dofs,1);
        L(esns.dofsg(end)) = -1/abs(sum(F));
        sy = zeros(str.dofs,1);
        sy(esns.fix_eq) = L(esns.fix_eq);
        sy(esns.free_eq) = -transpose(K(esns.free_eq,esns.free_eq))\(transpose(K(esns.fix_eq,esns.free_eq))*L(esns.fix_eq));
        elemK = elementK( str.L , str.H , str.Nx , str.Ny , str.T , str.E , str.nu );
        for e = 1:str.Nx*str.Ny
            eqn_num = esns.edof(e,:);
            g(e) = transpose(sy(eqn_num))*(p*rho(e)^(p-1) * elemK )*d(eqn_num);
        end
        
        %   Reverse filterin op
        g = w*g;
    end
end

%  Compliance constraint
function [c,ceq,gc,gceq] = constrFcn(rho,w,p,str,esns)
%     tic
    rho = w*rho;
    [d,K,F] = femsolve(rho,p,str,esns);
    F(esns.fix_eq) = K(esns.fix_eq,esns.free_eq)*d(esns.free_eq);
    c = F'*d - 5e-3;
%     toc
    if nargout > 1
        gc=zeros(str.Nx*str.Ny,1);
        elemK = elementK(str.L,str.H,str.Nx,str.Ny,str.T,str.E,str.nu);
        for e = 1:str.Nx*str.Ny
            eqn_num = list_dofs(e,str.Nx);
            gc(e) = -p*(rho(e))^(p-1)*(d(eqn_num)')*elemK*d(eqn_num);
        end
        gc = w*gc;
    end
    ceq=[];
    gceq=[];
end

%{
function [c,ceq,gc,gceq] = constrFcn(rho,w,p,str,esns)
    tic
    %--------------------------------------------------------------------------
    %   Apply density filter
    %--------------------------------------------------------------------------
    rho = w*rho;
    %--------------------------------------------------------------------------
    %   Solve FEM problem w filtered densities
    %--------------------------------------------------------------------------
    K = globalK(str.L,str.H,str.Nx,str.Ny,esns.edof,str.T,rho,p,str.E,str.nu);
    F = globalF_proj( str.L , str.H , str.Nx , str.Ny , str.T , str.P0 );
    d=sparse(str.dofs,1);
    d(esns.free_eq) = K(esns.free_eq,esns.free_eq)\F(esns.free_eq);
    %--------------------------------------------------------------------------
    %   Force an upper limit on the vertical disp at the input
    %--------------------------------------------------------------------------
    inputn = 1+str.Ny*(str.Nx+1);   %   Node number at the input (top left)
    c = -d(2*inputn)-0.01;     %   At most 0.5cm deflection
    toc
    %--------------------------------------------------------------------------
    %   Adjoint sensitivity for constraint grad
    %--------------------------------------------------------------------------
    if nargout > 1
        L = sparse(str.dofs,1);
        sy = zeros(str.dofs,1);
        L(2*inputn) = -1;
        
        sy(esns.free_eq) = -transpose(K(esns.free_eq,esns.free_eq))\L(esns.free_eq);
        elemK = elementK( str.L , str.H , str.Nx , str.Ny , str.T , str.E , str.nu );
        gc = zeros( str.Nx*str.Ny , 1 );
        for e = 1:str.Nx*str.Ny    
        
            eqn_num = list_dofs( e , str.Nx );
        
            gc(e) = transpose(sy(eqn_num)) * (p*rho(e)^(p-1) * elemK )  * d(eqn_num);
        
        end
    end

    %   Reverse filtering op
    gc = w*gc;
    ceq= [];
    gceq=[];
end
%}

function str = problem
    %   Bounding box dims (SI units)
    str.L = 0.16;   str.H = 0.04; str.T = 0.01;
    %   Gripper region dims (SI units)
    str.Lg = 0.025; str.Hg = 0.01;
    %   Material props
    str.E = 11e6;   str.nu = 0.4;
    %   Bounding box element numbers
    str.Nx = 128;    str.Ny = 64;
    %   total dofs
    str.dofs = 2*(str.Nx+1)*(str.Ny+1);
    %   Mesh sizes
    str.hx = str.L/str.Nx ; str.hy = str.H/str.Ny;
    %   Number of elements in gripping region
    str.Ngx = str.Lg/str.hx;        str.Ngy = str.Hg/str.hy;
    %   Uniform distributed force;
    str.P0 = -1e3;
end

function res = connect(str)
    %   Element dofs
    res.edof = list_dofs(1:str.Nx*str.Ny,str.Nx);

    %   Identify elements within gripping region
    bot = str.Nx:-1:str.Nx-str.Ngx+1;   %   Bottom row in gripper region
    res.Eg = [];                %   List elems in gripper region
    for j = 1:str.Ngy
        res.Eg = [ res.Eg bot+ (j-1)*str.Nx];
    end
    
    %   Identify dofs on gripper surf tip
    ng = str.Nx+1+str.Ngy*(str.Nx+1);
    res.dofsg = sort([ 2*ng-1 , 2*ng ]);
    
    %   Identify vertical dofs on bottom surf       
    res.botdofs = 2 * ( 1:str.Nx+1 );
    
    %   Segregate prescribed and free dofs
    res.fix_eq = sort([ res.botdofs , res.dofsg ]);
    res.free_eq = setdiff(1:2*(str.Nx+1)*(str.Ny+1),res.fix_eq); 
end
