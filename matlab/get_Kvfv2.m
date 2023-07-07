function [Kv, fv] = get_Kvfv2(tr, CONSTS)
    % Return contribution of second integral in assignment to K_v and f_v with linearized R_v
    M = size(tr.Points,1);
    num_triangles = size(tr.ConnectivityList,1);
    Kv = spalloc(M,M,7*M);
    fv = zeros(M,1);
    for k = 1:num_triangles
        inds = tr.ConnectivityList(k,:); % Indices of triangle points
        rz = tr.Points(inds,:);          % Triangle points
        for j = 1:3
            ind_j = inds(j);
            res_f = IntegralContainer.integral2d3(j,rz);
            fv(ind_j) = fv(ind_j) + CONSTS.V_mfv*res_f;
            for i = 1:3
                res = IntegralContainer.integral2d2(i,j,rz);
                ind_i = inds(i);
                Kv(ind_j, ind_i) = Kv(ind_j, ind_i) - CONSTS.r_q*CONSTS.V_mu/CONSTS.K_mu*res;
            end
        end
    end
end