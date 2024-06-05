function F = globalF_proj2( L , H , Nx , Ny , t , P0 )

    assert(isscalar(t), 'Uniform thickness needs to be assumed for force vector contributions')

    dofs = 2*(Nx+1)*(Ny+1);
    F = sparse(dofs,1);

    %   Iterate over elements on the LEFT FACE only!
    elems = 1:Nx:(1+(Ny-1)*Nx);
    for e = elems
        
        %   List global nodes for current elem
        %   Use only to size element force vector
        f_eq = list_dofs( e , Nx );
        
        %   Compute local F vector
        elemF = zeros(length(f_eq),1);
        detfactor2 = H/Ny*1/2;
        qp = 1/sqrt(3);
%         if rem(e,Nx)>0
%             E = mod(e,Nx);
%         else
%             E = Nx;
%         end
        E = find(elems==e);
        for yi = [qp -qp]
            h = H/Ny;
            ye = (E-1)*(1-yi)*h/2 + E*(1+yi)*h/2;
            qy = P0*(0.02 <= ye)*(ye <= 0.04);
            elemF(1) = elemF(1) + t*detfactor2*qy*(1-yi)/2;
            elemF(5) = elemF(5) + t*detfactor2*qy*(1+yi)/2;
        end

        F(f_eq) = F(f_eq) + elemF;
    end
end
