function K = get_K1(tr, sigma_r, sigma_z)
    % Return the contribution of the first integral in assignment to K_u or K_v
    M = size(tr.Points,1);
    num_triangles = size(tr.ConnectivityList,1);
    K = spalloc(M,M,7*M);
    for k = 1:num_triangles
        inds = tr.ConnectivityList(k,:); % Indices of triangle points
        rz = tr.Points(inds,:);          % Triangle points
        for i = 1:3
            for j = 1:3
                res = IntegralContainer.integral2d1(i,j,rz,sigma_r,sigma_z);
                ind_i = inds(i);
                ind_j = inds(j);
                K(ind_j, ind_i) = K(ind_j, ind_i) + res;
            end
        end
    end
end

