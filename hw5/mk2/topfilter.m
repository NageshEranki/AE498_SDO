function w = topfilter(L,H,Nx,Ny,ep)
    hx=L/Nx;hy=H/Ny;
    w=sparse(Nx*Ny,Nx*Ny);
    [xi,yi]=ndgrid(linspace(hx/2,L-hx/2,Nx),linspace(hy/2,H-hy/2,Ny));
    coords=[xi(:) yi(:)];
    for i=1:Nx*Ny
        for j=i:Nx*Ny
            w(i,j)=max(0,ep-norm(coords(i,:)-coords(j,:)));
        end
    end
    d=diag(w);
    w=w-diag(d);
    w=w+w'+diag(d);
    for i = 1:Nx*Ny
        w(i,:)=w(i,:)/sum(w(i,:));
    end
end
%{
clear;clc;
%--------------------------------------------------------------------------
%   Mesh set-up
%--------------------------------------------------------------------------
L=6;   H=2;
Nx=60;  Ny=20;
hx=L/Nx;hy=H/Ny;
%--------------------------------------------------------------------------
%   Filter radius
%--------------------------------------------------------------------------
ep = 0.15;
%--------------------------------------------------------------------------
rho=kron([ones(Ny/2,1); zeros(Ny/2,1)],ones(1,Nx));
% rho=kron(ones(Ny,1),1:Nx);
I=mat2gray(rho);
figure
imshow(I,'initialmagnification',1600)
[xi,yi]=ndgrid(linspace(hx/2,L-hx/2,Nx),linspace(hy/2,H-hy/2,Ny));
coords=[xi(:) yi(:)];
w=sparse(Nx*Ny,Nx*Ny);
tic
for i=1:Nx*Ny
    for j=i:Nx*Ny
        w(i,j)=max(0,ep-norm(coords(i,:)-coords(j,:)));
    end
end
toc
d=diag(w);
w=w-diag(d);
w=w+w'+diag(d);
for i = 1:Nx*Ny
    w(i,:)=w(i,:)/sum(w(i,:));
end
rho=reshape(transpose(flipud(rho)),Nx*Ny,1);
rho=w*rho;
rho=flipud(transpose(reshape(rho,Nx,Ny)));

I=mat2gray(rho);
figure
imshow(I,'initialmagnification',1600)
%}