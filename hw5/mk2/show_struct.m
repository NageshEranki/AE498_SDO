function show_struct(Nx,Ny,rho)
    rho=flipud(transpose(reshape(rho,Nx,Ny)));
    I=mat2gray(rho);
    figure
    imshow(I,'initialmagnification',1600)
%     colormap(gray);
%     imagesc(-rho);
%     axis equal
end
