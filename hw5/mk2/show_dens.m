function show_dens(rho,L,H,nx,ny)
    figure

    [xi,yi] = ndgrid( linspace(0,L,nx+1) , linspace(0,H,ny+1) );

    nodeMap = list_nodes(1:nx*ny,nx);

    X = zeros(4,nx*ny); Y = zeros(4,nx*ny);

    for e = 1:nx*ny
        X(:,e) = [ xi(nodeMap(e,1)) ; xi(nodeMap(e,2)) ; xi(nodeMap(e,4)) ; xi(nodeMap(e,3)) ];
        Y(:,e) = [ yi(nodeMap(e,1)) ; yi(nodeMap(e,2)) ; yi(nodeMap(e,4)) ; yi(nodeMap(e,3)) ];
    end

    patch(X,Y,rho,'EdgeColor','none')
    colormap('gray')
    axis equal
end