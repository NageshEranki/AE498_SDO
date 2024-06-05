function res = list_nodes( e , Nx )

    res = zeros(length(e),4);

    res(:,1) = ( mod(e,Nx) > 0 ).*( mod(e,Nx) + floor(e/Nx) * (Nx+1) ) + ...
               ( mod(e,Nx) == 0 ).*( floor(e/Nx) * (Nx+1)-1 );
    
    res(:,2) = res(:,1) + 1;
    
    res(:,3) = res(:,1) + Nx + 1;
    
    res(:,4) = res(:,2) + Nx + 1;

end