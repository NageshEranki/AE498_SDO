function F = globalF( L , H , Nx , Ny , t , P0 )

    assert(isscalar(t), 'Uniform thickness needs to be assumed for force vector contributions')

    dofs = 2*(Nx+1)*(Ny+1);
    %F = zeros(dofs,1);
    F = sparse(dofs,1);

    %   Iterate over elements on the top row only!
    for e = (Nx*(Ny-1)+1):(Nx*Ny)

        f_eq = list_dofs( e , Nx );
        
        Felem = elementF( L , H , Nx , Ny , t , P0 , e );

        F(f_eq) = F(f_eq) + Felem;
    end
end
