%function Ke = elementK(L , H , Nx , Ny , t , r , p , E , nu)
function Ke = elementK(L , H , Nx , Ny , t , E , nu)
    
    %assert(isscalar(r) , "Element density must be scalar for elementK()")
    
    Ke = zeros(8,8);

    %   Constitutive matrix for plane stress quad element
    %   Using SIMP penalization
%     p = 3;
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];

    %   Quadrature point
    qp = 1/sqrt(3);

    %   Inverse of the jacobian;    will be a constant matrix
    jacinv = [2*Nx/L 0; 0 2*Ny/H];

    %   Integration transformation rule
    detfactor = L/Nx * H/Ny * 1/4;

    for xi = [-qp qp]

        for eta = [-qp qp]

        %   Shape function gradient on ref element
        lgrad = [-(1-eta)/4 (1-eta)/4 -(1+eta)/4 (1+eta)/4; -(1-xi)/4 -(1+xi)/4 (1-xi)/4 (1+xi)/4];

        
        %   Tranform local gradients
        ggrad = zeros(size(lgrad));

        for j = 1:size(lgrad,2)
            ggrad(:,j) = jacinv*lgrad(:,j);
        end


        %   Assemble B; strain interpolation matrix
        B = [];
        for j = 1:size(lgrad,2)
            B = [B [ggrad(1,j) 0; 0 ggrad(2,j); ggrad(2,j) ggrad(1,j)]];
        end

        
        Ke = Ke + t*detfactor*B'*C*B;

        end
    end

end
