function [K, f] = get_Kf3(tr, rho, C_amb)
    % Return the contribution of the third integral in assignment to K_u/K_v and f_u/f_v
    M = size(tr.Points,1);
    % Determine triangle edges on Gamma_2 (outer boundary)
    fb = freeBoundary(tr);
    r = tr.Points(:,1);
    r_vals_edges = r(fb);
    r_equal_zero = r_vals_edges==0; % For
    edge_on_outer_edge = ~(r_equal_zero(:,1) & r_equal_zero(:,2));
    outer_edges = fb(edge_on_outer_edge,:); % [(ind1,ind2)] of edges on Gamma_2

    num_edges = size(outer_edges,1);
    K = spalloc(M,M,7*M);
    f = zeros(M,1);
    for k = 1:num_edges
        inds = outer_edges(k,:); % Indices of edge extremities
        rz = tr.Points(inds,:);  % Corresponding points
        for j = 1:2
            ind_j = inds(j);
            f(ind_j) = f(ind_j) + rho*C_amb*IntegralContainer.integral1d1(j,rz);
            for i = 1:2
                res = IntegralContainer.integral1d2(i,j,rz);
                ind_i = inds(i);
                K(ind_j,ind_i) = K(ind_j,ind_i) + rho*res;
            end
        end
    end
end