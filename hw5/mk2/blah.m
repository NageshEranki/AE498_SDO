function res = blah(Nx,Ny)
    L = 16; H = 4;
    l = 2.5; h = 1;
    
    hx = L/Nx;  hy = H/Ny;
    res = [];
    for i = 1:h/hy
        res = [res i*Nx:-1:(i*Nx-l/hx)+1];
    end
    
end

