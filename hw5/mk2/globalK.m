function K = globalK( L , H , Nx , Ny , edof , t , r , p , E , nu )

	assert( ( length(r) == Nx * Ny ) , "List of individual densities is incomplete")

	dofs = 2 * ( Nx + 1 ) * ( Ny + 1 );

	%K = zeros(dofs,dofs);
	K = sparse(dofs,dofs);

	%	Compute local K once;
	%	Impose SIMP penalty during assembly
	Kelem = elementK( L , H , Nx , Ny , t , E , nu );

% 	edof = list_dofs(1:Nx*Ny,Nx);

	for e = 1:Nx*Ny

		%       Stitch local K into global K
		K(edof(e,:),edof(e,:)) = K( edof(e,:),edof(e,:) ) + (r(e)^p)*Kelem;

	end

end
