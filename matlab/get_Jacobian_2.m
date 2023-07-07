function [J, h] = get_Jacobian_2(order, tr, CONSTS, c)
% Returns the linearization of the nonlinear function H with the given order,
% given a point c.
% Returns the jacobian J and the evaluation of the function H in point c.
    M = size(tr.Points,1);
    num_triangles = size(tr.ConnectivityList,1);
    J1 = spalloc(M,M,7*M);
    J2 = spalloc(M,M,7*M);
    J3 = spalloc(M,M,7*M);
    J4 = spalloc(M,M,7*M);
    hu = zeros(M, 1); 
    hv = zeros(M, 1);
    
    [quadrature_points, w] = get_quadrature_points_weights(order);
    
    N = zeros(length(w),3);
    for k = 1:length(w)
        % Each row contains the basis functions evaluated at quadrature point
        N(k,:) = [N1(quadrature_points(k,:)) N2(quadrature_points(k,:)) N3(quadrature_points(k,:))];
    end

    for k = 1:num_triangles
        inds = tr.ConnectivityList(k,:); % Indices of triangle points
        rz = tr.Points(inds,:);                     % Triangle points
        r = rz(:,1);
        det_J = det_jacobian(rz);
        cu = c(inds); %% Coefficients with active basis functions of C_u
        cv = c(M + inds); %%Coefficients with active basis functions of C_v

        Nr = N*r;
        Ncu = N*cu;
        Ncv = N*cv;
        
        for m=1:3
            ind_m = inds(m);
            Nm = N(:, m); % Basis function m evaluated at quadrature points.
            
            noemer1 = 1+1/CONSTS.K_mv.*Ncv; %First term in denominator of alpha
            noemer2 = CONSTS.K_mu + Ncu; %Second term in denominator of alpha
            noemer3 = (1 + 1/CONSTS.K_mfu.*Ncu).^2; % Denominator of gamma

            delta_integrand = (Nr).*Ncu.*Nm./(noemer1)./(noemer2);
            delta = CONSTS.V_mu*det_J*(w'*delta_integrand);

            epsilon_integrand = (Nr).*Nm./(1+1/CONSTS.K_mfu.*Ncu);
            epsilon = CONSTS.V_mfv*det_J*(w'*epsilon_integrand);

            hu(ind_m) = hu(ind_m) + delta;
            hv(ind_m) = hv(ind_m) - CONSTS.r_q * delta - epsilon;

            for i=m:3
                ind_i = inds(i);
                Ni = N(:, i); % Basis function i evaluated at quadrature points.
                
                tripleN = Nr.*Ni.*Nm;

                alpha_integrand = tripleN./noemer1./(noemer2.^2);
                alpha = CONSTS.V_mu*CONSTS.K_mu*det_J*(w'*alpha_integrand);

                beta_integrand = Ncu.*tripleN./(noemer1.^2)./noemer2;
                beta = -CONSTS.V_mu/CONSTS.K_mv*det_J*(w'*beta_integrand);

                gamma_integrand = tripleN./noemer3;
                gamma = CONSTS.V_mfv/CONSTS.K_mfu*det_J*(w'*gamma_integrand);

                J1(ind_m, ind_i) = J1(ind_m, ind_i) + alpha;
                J2(ind_m, ind_i) = J2(ind_m, ind_i) + beta;
                J3(ind_m, ind_i) = J3(ind_m, ind_i) - CONSTS.r_q*alpha + gamma;
                J4(ind_m, ind_i) = J4(ind_m, ind_i) - CONSTS.r_q*beta;

                if m ~= i %% Exploit symmetry
                    J1(ind_i, ind_m) = J1(ind_i, ind_m) + alpha;
                    J2(ind_i, ind_m) = J2(ind_i, ind_m) + beta;
                    J3(ind_i, ind_m) = J3(ind_i, ind_m) - CONSTS.r_q*alpha + gamma;
                    J4(ind_i, ind_m) = J4(ind_i, ind_m) - CONSTS.r_q*beta;
                end
            end
        end
    end

J = [J1, J2; J3, J4];
h = [hu; hv];
end

function J_det = det_jacobian(rz)
    r = rz(:,1); z = rz(:,2);
    J = [r(2)-r(1) r(3)-r(1); z(2)-z(1) z(3)-z(1)];
    J_det = abs(det(J));
end
function res = N1(xi_eta)
    res = 1-xi_eta(1)-xi_eta(2);
end
function res = N2(xi_eta)
    res = xi_eta(1);
end
function res = N3(xi_eta)
    res = xi_eta(2);
end
