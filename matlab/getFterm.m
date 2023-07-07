function f = getFterm(tr, consts, ms)

    Ru = @(r, z) consts.V_mu*ms.Cu(r, z)./(consts.K_mu + ms.Cu(r,z))./(1. + ms.Cv(r,z)./consts.K_mv);
    Rv = @(r, z) consts.r_q*Ru(r,z) + consts.V_mfv./(1. + ms.Cu(r,z)./consts.K_mfu);
    
    Qu = @(r, z) consts.sigma_ur*(ms.dCu_dr(r,z) + r.*ms.dCu_dr2(r,z)) ...
               + consts.sigma_uz*(r.*ms.dCu_dz2(r,z)) ...
               - r.*Ru(r, z);
    
    Qv = @(r, z) consts.sigma_vr*(ms.dCv_dr(r,z) + r.*ms.dCv_dr2(r,z)) ...
               + consts.sigma_vz*(r.*ms.dCv_dz2(r,z)) ...
               + r.*Rv(r, z);
    
    nb_of_triangles = size(tr.ConnectivityList, 1);
    M = size(tr.Points, 1);

    f = zeros(2*M, 1);

    N1 = @(xi, eta) 1-xi-eta;
    N2 = @(xi, eta) xi;
    N3 = @(xi, eta) eta;
    NN = {N1, N2, N3};
    
    %% Production term contribution
    for k=1:nb_of_triangles
        inds = tr.ConnectivityList(k, :);
        rz = tr.Points(inds, :);
        r = rz(:,1); z = rz(:,2);
        J = [r(2)-r(1) r(3)-r(1); z(2)-z(1) z(3)-z(1)];
        det_j = abs(det(J));
        for m=1:3
            ind_m = inds(m);
            Nm = NN{m};
            rf = @(xi, eta) r(1)*N1(xi, eta) + r(2)*N2(xi, eta) + r(3)*N3(xi, eta);
            zf = @(xi, eta) z(1)*N1(xi, eta) + z(2)*N2(xi, eta) + z(3)*N3(xi, eta);
            int1 = @(xi, eta) Qu(rf(xi, eta), zf(xi, eta)).*Nm(xi, eta);
            int2 = @(xi, eta) Qv(rf(xi, eta), zf(xi, eta)).*Nm(xi, eta);
            
            eta_max = @(xi) 1-xi;
            f(ind_m) = f(ind_m) - det_j*integral2(int1, 0, 1, 0, eta_max);
            f(ind_m + M) = f(ind_m + M) - det_j*integral2(int2, 0, 1, 0, eta_max);
        end
    end
    
    %%Boundery term contribution
    outer_edges = outer_boundary(tr);

    num_edges = size(outer_edges,1);
    for k = 1:num_edges
        inds = outer_edges(k,:); % Indices of edge extremities
        rz = tr.Points(inds,:);  % Corresponding points
        assert(all(size(rz) == [2, 2])); 
        r = rz(:,1); z = rz(:,2);
        deltar = r(2) - r(1);
        deltaz = z(2) - z(1);
        E2 = sqrt(deltar^2 + deltaz^2);
        normalr = deltaz/E2;
        normalz = -deltar/E2;
        assert(abs((deltar*normalr + deltaz*normalz)) <= 1e-14);
        assert(abs(norm([normalr, normalz])) - 1 <= 1e-14);

        for j = 1:2
            ind_j = inds(j);
            Nj = NN{j};
            rf = @(xi) r(1)*N1(xi, 0) + r(2)*N2(xi, 0);
            zf = @(xi) z(1)*N1(xi, 0) + z(2)*N2(xi, 0);
            int1 = @(xi) Nj(xi, 0).*( ...
                (normalr.*rf(xi).*consts.sigma_ur.*ms.dCu_dr(rf(xi), zf(xi))) ...
              + (normalz.*rf(xi).*consts.sigma_uz.*ms.dCu_dz(rf(xi), zf(xi))) ...
             );
            int2 = @(xi) Nj(xi, 0).*( ...
                (normalr.*rf(xi).*consts.sigma_vr.*ms.dCv_dr(rf(xi), zf(xi))) ...
              + (normalz.*rf(xi).*consts.sigma_vz.*ms.dCv_dz(rf(xi), zf(xi))) ...
             );
            f(ind_j) = f(ind_j) + E2*integral(int1, 0, 1);
            f(ind_j + M) = f(ind_j + M) + E2*integral(int2, 0, 1);
        end
    end
end
