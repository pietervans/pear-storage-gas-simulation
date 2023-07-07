function h = compute_h(order, tr, CONSTS, c)
    M = size(tr.Points,1);
    num_triangles = size(tr.ConnectivityList,1);
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

            delta_integrand = (Nr).*Ncu.*Nm./(noemer1)./(noemer2);
            delta = CONSTS.V_mu*det_J*(w'*delta_integrand);

            epsilon_integrand = (Nr).*Nm./(1+1/CONSTS.K_mfu.*Ncu);
            epsilon = CONSTS.V_mfv*det_J*(w'*epsilon_integrand);

            hu(ind_m) = hu(ind_m) + delta;
            hv(ind_m) = hv(ind_m) - CONSTS.r_q * delta - epsilon;
        end
    end
    h = [hu;hv];
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
