clear;clc

str = problem;
esns = connect(str);
globK = @(rho) globalK(str.L,str.H,str.Nx,str.Ny,esns.edof,str.T,rho,3,str.E,str.nu);
rho = 0.1*ones(str.Nx*str.Ny,1);
rho(esns.Eg) = 1e-3;

p=3;
fdfunc = @(rho) findiff(rho,p);
adjfunc= @(rho) adjmethod(rho,p);

s1=fdfunc(rho);s2=adjfunc(rho);
max(abs(s1-s2))

function dfdrho = adjmethod(rho,p)
    str=problem;
    esns=connect(str);
    
    dfdrho = zeros(str.Nx*str.Ny,1);
    
    K = globalK(str.L,str.H,str.Nx,str.Ny,esns.edof,...
                str.T, rho, p , str.E, str.nu);
    F = globalF_proj(str.L,str.H,str.Nx,str.Ny,str.T,str.P0);
    d = sparse(str.dofs,1);
    d(esns.free_eq) = K(esns.free_eq,esns.free_eq)\F(esns.free_eq);  
    F(esns.fix_eq) = K(esns.fix_eq,esns.free_eq)*d(esns.free_eq);
    
    L = sparse(str.dofs,1);
    L(esns.dofsg(end)) = 1;
    sy = zeros(str.dofs,1);
    sy(esns.fix_eq) = L(esns.fix_eq);
    sy(esns.free_eq) = -transpose(K(esns.free_eq,esns.free_eq))\(transpose(K(esns.fix_eq,esns.free_eq))*L(esns.fix_eq));
    elemK = elementK( str.L , str.H , str.Nx , str.Ny , str.T , str.E , str.nu );
    for e = 1:str.Nx*str.Ny
%         if ~ismember(e,esns.Eg)
        if 1
            eqn_num = esns.edof(e,:);
            dfdrho(e) = transpose(sy(eqn_num))*(p*rho(e)^(p-1) * elemK )*d(eqn_num);
        end
    end
end

function dfdrho = findiff(rho,p)
    str = problem;
    esns = connect(str);
    h = 1e-10;
    
    K0 = globalK(str.L,str.H,str.Nx,str.Ny,esns.edof,...
                str.T, rho, p, str.E, str.nu);
    F0 = globalF_proj(str.L,str.H,str.Nx,str.Ny,str.T,str.P0);
    d0 = sparse(str.dofs,1);
    d0(esns.free_eq) = K0(esns.free_eq,esns.free_eq)\F0(esns.free_eq);
    F0(esns.fix_eq) = K0(esns.fix_eq,esns.free_eq)*d0(esns.free_eq);

    f0 = F0(esns.dofsg(end));
    
    dfdrho = zeros(str.Nx*str.Ny,1);
    F = globalF_proj(str.L,str.H,str.Nx,str.Ny,str.T,str.P0);
    for e = 1:str.Nx*str.Ny
%         if ~ismember(e,esns.Eg) 
        if 1
            rho2 = rho;
            rho2(e) = rho2(e) + h;
            K = globalK(str.L,str.H,str.Nx,str.Ny,esns.edof,...
                    str.T, rho2, p , str.E, str.nu);
            d = sparse(str.dofs,1);
            d(esns.free_eq) = K(esns.free_eq,esns.free_eq)\F(esns.free_eq);
            Fi = F;
            Fi(esns.fix_eq) = K(esns.fix_eq,esns.free_eq)*d(esns.free_eq);

            fi = Fi(esns.dofsg(end));

            dfdrho(e) = (fi-f0)/h;
        end
    end
end

function str = problem
    %   Bounding box dims (SI units)
    str.L = 0.16;   str.H = 0.04; str.T = 0.01;
    %   Gripper region dims (SI units)
    str.Lg = 0.025; str.Hg = 0.01;
    %   Material props
    str.E = 11e6;   str.nu = 0.4;
    %   Bounding box element numbers
    str.Nx = 32;    str.Ny = 8;
    %   total dofs
    str.dofs = 2*(str.Nx+1)*(str.Ny+1);
    %   Mesh sizes
    str.hx = str.L/str.Nx;  str.hy = str.H/str.Ny;
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