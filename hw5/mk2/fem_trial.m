%--------------------------------------------------------------------------
% Solves the FEM problem for HW 5 set-up
%		
%	MAKING USE OF THE SYMMETRY OF THE PROBLEM
%
%--------------------------------------------------------------------------
clear;clc
%--------------------------------------------------------------------------
%   Material constants (in SI units)
%--------------------------------------------------------------------------
E = 75e6;   nu = 0.3;   t = 0.1; p = 3;
P0 = -100e3;
%--------------------------------------------------------------------------
%   Mesh params
%--------------------------------------------------------------------------
xmin = 0;   xmax = 6;  ymin = 0;   ymax = 2;
Nx = 60;   Ny = 20;
%--------------------------------------------------------------------------
%   Time start
tic
%--------------------------------------------------------------------------
%   Total number of degrees of freedom
dofs = 2*(Nx+1)*(Ny+1);
%   Elements densities
r = 0.4*ones(Nx*Ny,1);
edof = list_dofs(1:Nx*Ny,Nx);
%--------------------------------------------------------------------------
%   Assemble global K
%--------------------------------------------------------------------------
K = globalK( xmax-xmin , ymax-ymin , Nx , Ny , edof ,  t ,  r , p , E , nu );
%--------------------------------------------------------------------------
%   Global force vector
%--------------------------------------------------------------------------
F = globalF( xmax-xmin , ymax-ymin , Nx , Ny , t , P0 );
%--------------------------------------------------------------------------
%	Imposing prescribed disp.
%	
%	Left surface:	Fix support
%	Right surface:	Roller support
%--------------------------------------------------------------------------
leftsurf = 1:Nx+1:(Nx+1)*(Ny+1);
rightsurf = Nx+1:Nx+1:(Nx+1)*(Ny+1);
fix_eq = sort( [ 2*leftsurf-1 2*leftsurf 2*rightsurf-1 ] );
%--------------------------------------------------------------------------
%   List eqn numbers corresp. to free dofs
free_eq=setdiff(1:dofs,fix_eq);

d = zeros(dofs,1);
d(free_eq) = K(free_eq,free_eq)\F(free_eq);

toc

strainE_density = 0.5*d'*K*d/(xmax*ymax*t)

%   Plot results
[xi,yi] = ndgrid(linspace(xmin,xmax,Nx+1),linspace(ymin,ymax,Ny+1));
xi = xi(:); yi = yi(:);
d2 = (0.1*(xmax-xmin))/max(abs(d),[],'all')*d;
%d2 = d;

dx = d2(1:2:end);    dy = d2(2:2:end);
plot(xi+dx,yi+dy,'xb')
axis equal
% saveas(gcf,'~/femv2.png')
