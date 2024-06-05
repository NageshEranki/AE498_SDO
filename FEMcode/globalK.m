function K = globalK( L , H , Nx , Ny , t , E , nu )

    if( ~isscalar(t) )

        assert( ( length(t) == Nx * Ny ) , "List of individual thicknesses is incomplete")

    end

    dofs = 2 * ( Nx + 1 ) * ( Ny + 1 );

    K = zeros(dofs,dofs);

    for e = 1:Nx*Ny

        %       Compute local K each time; t may vary for each element
        if(isscalar(t))

            Kelem = elementK( L , H , Nx , Ny , t , E , nu );

        else

            Kelem = elementK( L , H , Nx , Ny , t(e) , E , nu );

        end

        %       List the global dofs for the current elem
        K_eq = list_dofs( e , Nx );
        
        %       Stitch local K into global K
        K(K_eq,K_eq) = K( K_eq , K_eq ) + Kelem;

    end

end