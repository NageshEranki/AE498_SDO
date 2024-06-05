% Function generates a coefficient matrix used for filtering of element
% densities
% W contains the weight coefficients for a linear filter with radius given
% by the parameter 'rad'
% Inputs:
% Ex = number of elements in the x-direction
% Ey = number of elements in the y-direction
% nh = neighborhood size (approximate radius neighborhood measured in elements)
% X  = x-coordinates of all nodes in mesh
% Y  = y-coordinates of all nodes in mesh
% nodeMap = array indicating which nodes belong to each element

function W = density_filter(Ex, Ey, nh, X, Y, nodeMap)

% compute element centroid locations

n_e = Ex*Ey; % Number of elements in the mesh
X_cent = zeros(n_e, 1);
Y_cent = zeros(n_e, 1);

for i = 1:n_e
    
    X_cent(i) = sum(X(nodeMap(i, :)))/4;
    Y_cent(i) = sum(Y(nodeMap(i, :)))/4;
    
end    

% find 'neighborhood' size (i.e. radius (measured in elements) of the filter neighborhood)
% based on characteristic length of elements
scale = max(X_cent(2)-X_cent(1), Y_cent(2)-Y_cent(1));
rad = nh*scale;

W = zeros(n_e, n_e);

for i = 1:n_e

    weight_sum = 0;

    for iw = 1:n_e

        dist = sqrt((X_cent(i) - X_cent(iw))^2 + (Y_cent(i) - Y_cent(iw))^2);
        coeff = rad - dist;

        if coeff > 0
            W(i, iw) = coeff;
            weight_sum = weight_sum + coeff;
        end

    end

    % normalize weights
    W(i, :) = W(i, :)/weight_sum;

end