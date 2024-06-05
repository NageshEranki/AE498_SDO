function elemF = elementF( L , H , Nx , Ny , t , P0 , e )
    
    %   List global nodes for current elem
    f_eq = list_dofs( e , Nx );
    
    %   Compute local F vector
    elemF = zeros(length(f_eq),1);
    detfactor2 = L/Nx*1/2;
    qp = 1/sqrt(3);
    for xi = [qp -qp]
        h = L/Nx;
        xe = (e-1)*(1-xi)*h/2 + e*(1+xi)*h/2;
        qy = P0*(L-0.1 <= xe)*(xe <= L);
        elemF(2) = elemF(2) + t*detfactor2*qy*(1-xi)/2;
        elemF(4) = elemF(4) + t*detfactor2*qy*(1+xi)/2;
    end
    
end