function Ku = get_Ku2(tr, CONSTS)
    % Return contribution of second integral in assignment to K_u with linearized R_u
    M = size(tr.Points,1);
    num_triangles = size(tr.ConnectivityList,1);
    Ku = spalloc(M,M,7*M);
    for k = 1:num_triangles
        inds = tr.ConnectivityList(k,:); % Indices of triangle points
        rz = tr.Points(inds,:);          % Triangle points
        for i = 1:3
            for j = 1:3
                res = IntegralContainer.integral2d2(i,j,rz);
                ind_i = inds(i);
                ind_j = inds(j);
                Ku(ind_j, ind_i) = Ku(ind_j, ind_i) + CONSTS.V_mu/CONSTS.K_mu*res;
            end
        end
    end
end