function [ d , K , F , fix_eq , free_eq ] = fem_pstress( L , H  , Nx , Ny , t , E , nu , P0 )
    
    T = 0.05;

%     sprintf('Using a thickness of %f for force vector calcs',T)

    K = globalK( L , H , Nx , Ny , t , E , nu );

    F = globalF( L , H , Nx , Ny , T , P0 );

   
    fix_n = 1 + (Nx+1)*( (1:Ny+1) - 1);
    fix_eq = zeros(1,2*length(fix_n));
    for i = 1:length(fix_n)
        fix_eq(2*i-1) = 2*fix_n(i)-1; fix_eq(2*i) = 2*fix_n(i);
    end

    %   List eqn numbers corresp. to free dofs
    dofs = 2*(Nx+1)*(Ny+1);
    free_eq=setdiff(1:dofs,fix_eq);
    
    %   extract rows and columns of the free dofs
    Kred=K(free_eq,free_eq);
    Fred=F(free_eq);
    
    dred=Kred\Fred;
    
    
    %   Rebuild d
    d=zeros(dofs,1);
    for i=1:length(free_eq)
        d(free_eq(i))=dred(i);
    end
    
end