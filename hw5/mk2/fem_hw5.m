function [ d , K , F , fix_eq , free_eq ] = fem_hw5( L , H  , Nx , Ny , edof , t , r , p , E , nu , P0 )
	
	%----------------------------------------------------------------------
	%	Solves the FEA model for hw5
	%	
	%	MAKING USE OF SYMMETRY OF THE PROBLEM
	%----------------------------------------------------------------------
	%   Total number of degrees of freedom
	dofs = 2*(Nx+1)*(Ny+1);
	%----------------------------------------------------------------------
	%   Assemble global K
	%----------------------------------------------------------------------
	K = globalK( L , H , Nx , Ny , edof , t ,  r , p , E , nu );
	%----------------------------------------------------------------------
	%   Global force vector
	%----------------------------------------------------------------------
	F = globalF( L , H , Nx , Ny , t , P0 );
	%--------------------------------------------------------------------------
	%	Imposing prescribed disp.
	%	
	%	Left surface:	Fix support; no x,y disp.
	%	Right surface:	Roller support; no x disp.
	%--------------------------------------------------------------------------
	leftsurf = 1:Nx+1:(Nx+1)*(Ny+1);
	rightsurf = Nx+1:Nx+1:(Nx+1)*(Ny+1);
	fix_eq = sort( [ 2*leftsurf-1   2*leftsurf   2*rightsurf-1 ] );
	%----------------------------------------------------------------------
	%   List eqn numbers corresp. to free dofs
	free_eq=setdiff(1:dofs,fix_eq);

	%   Only free dofs have non-zero disp.
	d = zeros(dofs,1);
	d(free_eq) = K(free_eq,free_eq)\F(free_eq);

end
