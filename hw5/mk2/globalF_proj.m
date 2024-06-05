function F = globalF_proj( L , H , Nx , Ny , t , P0 )

    assert(isscalar(t), 'Uniform thickness needs to be assumed for force vector contributions')

    dofs = 2*(Nx+1)*(Ny+1);
    F = sparse(dofs,1);

    %   Iterate over elements on the top row only!
    for e = (Nx*(Ny-1)+1):(Nx*Ny)

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
            qy = P0*(0 <= xe)*(xe <= 0.04);
            elemF(6) = elemF(6) + t*detfactor2*qy*(1-xi)/2;
            elemF(8) = elemF(8) + t*detfactor2*qy*(1+xi)/2;
        end

        F(f_eq) = F(f_eq) + elemF;
    end
end
