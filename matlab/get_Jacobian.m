function [J, h] = get_Jacobian(tr, CONSTS, c)
%Returns the linearization of the nonlinear function H,
% given a point c.
% Returns the jacobian J and the evaluation of the function H in point c.
    M = size(tr.Points,1);
    num_triangles = size(tr.ConnectivityList,1);
    J1 = zeros(M);
    J2 = zeros(M);
    J3 = zeros(M);
    J4 = zeros(M);
    hu = zeros(M, 1); 
    hv = zeros(M, 1);

    for k = 1:num_triangles
        inds = tr.ConnectivityList(k,:); % Indices of triangle points
        for m = 1:3
            ind_m = inds(m);
            %% Contributions to Hu(ck) and Hv(ck)
            delta = IntegralContainer.compute_delta(m,tr,k,c,CONSTS);
            epsilon = IntegralContainer.compute_epsilon(m, tr, k, c, CONSTS);
            hu(ind_m) = hu(ind_m) + delta;
            hv(ind_m) = hv(ind_m) - CONSTS.r_q * delta - epsilon;
            for i = 1:3
                ind_i = inds(i);
                alpha = IntegralContainer.compute_alpha(i, m, tr, k, c, CONSTS);
                beta = IntegralContainer.compute_beta(i, m, tr, k, c, CONSTS);
                gamma = IntegralContainer.compute_gamma(i, m, tr, k, c, CONSTS);

                J1(ind_m, ind_i) = J1(ind_m, ind_i) + alpha;
                J2(ind_m, ind_i) = J2(ind_m, ind_i) + beta;
                J3(ind_m, ind_i) = J3(ind_m, ind_i) - CONSTS.r_q*alpha + gamma;
                J4(ind_m, ind_i) = J4(ind_m, ind_i) - CONSTS.r_q*beta;
            end
        end
    end

J = [J1, J2; J3, J4];
h = [hu; hv];
end