  function elemF = elementF( L , H , Nx , Ny , t , P0 , e )
    
    %   List global nodes for current elem
    %   Use only to size element force vector
    f_eq = list_dofs( e , Nx );
    
    %   Compute local F vector
    elemF = zeros(length(f_eq),1);
    detfactor2 = L/Nx*1/2;
    qp = 1/sqrt(3);
    if rem(e,Nx)>0
        E = mod(e,Nx);
    else
        E = Nx;
    end
    for xi = [qp -qp]
        h = L/Nx;
        xe = (E-1)*(1-xi)*h/2 + E*(1+xi)*h/2;
        qy = P0*(L-0.25 <= xe)*(xe <= L);
        elemF(6) = elemF(6) + t*detfactor2*qy*(1-xi)/2;
        elemF(8) = elemF(8) + t*detfactor2*qy*(1+xi)/2;
    end
    
end
