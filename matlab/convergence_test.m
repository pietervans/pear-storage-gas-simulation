%% Convergence test
CONSTS = get_constants('orchard');
element_sizes = [5,3,2,1,0.5];
mesh_list = ["uniform_5mm";
             "uniform_3mm";
             "uniform_2mm";
             "uniform_1mm";
             "uniform_0p5mm"];
quadrature_order = 3;

% %%Chosen solution
ms = struct();
ms.Cu = @(r, z) 0.05 + 0.05*sin(10*r) + z.^2;
ms.Cv = @(r, z) 2 + 3*r.^2 + sin(15*z);
    
ms.dCu_dr = @(r, z) 0.05*10*cos(10*r);
ms.dCu_dz = @(r, z) 2*z;
ms.dCv_dr = @(r, z) 6*r;
ms.dCv_dz = @(r, z) 15*cos(15*z);

ms.dCu_dr2 = @(r, z) -10*10*0.05*sin(10*r);
ms.dCu_dz2 = @(r, z) 2;
ms.dCv_dr2 = @(r, z) 6;
ms.dCv_dz2 = @(r, z) -15*15*sin(15*z);

results = zeros(size(mesh_list, 1), 2);
for mesh_i=1:size(mesh_list, 1)
    "Computing for " + mesh_list(mesh_i)

    tr = get_triangulation(mesh_list(mesh_i));
    
    tic()
    M = size(tr.Points,1);
    Ku = spalloc(M,M,7*M);
    Kv = spalloc(M,M,7*M);

    % Integral 1 in assignment (Equations (5) and (6))
    Ku = Ku + get_K1(tr, CONSTS.sigma_ur, CONSTS.sigma_uz);
    Kv = Kv + get_K1(tr, CONSTS.sigma_vr, CONSTS.sigma_vz);
    
    % Integral 3 in assigment, and production term, given manufactured solution
    f = getFterm(tr, CONSTS, ms);

    % Exact solution
    c_exact = zeros(2*M, 1);
    for i=1:M
        point = tr.Points(i, :);
        c_exact(i) = ms.Cu(point(1), point(2));
        c_exact(i+M) = ms.Cv(point(1), point(2));
    end

    c = getInitialC(Ku, Kv, f(1:M), f(M+1:end), tr, CONSTS, false);
    "Get initial c: " + toc()
    
    K = [Ku, sparse(M,M); sparse(M,M) Kv];
    h_c0 = compute_h(quadrature_order, tr, CONSTS, c);
    backwarderror = norm(K*c - f + h_c0)/norm(h_c0);
    
    forwardcu = norm(c(1:M)-c_exact(1:M))/norm(c_exact(1:M));
    forwardcv = norm(c(M:end)-c_exact(M:end))/norm(c_exact(M:end));
    
    it = 0;
    while backwarderror > 1e-13 && it < 10
        
        % Linearization of nonlinear function
        tic();
        [J, h_ck] = get_Jacobian_2(quadrature_order, tr, CONSTS, c);
        "Get Jacobian: " + toc();
        
        % Solve system
        A = K + J;
        b = f - h_ck + J*c;
        
        c = A\b;
       
        h = compute_h(quadrature_order, tr, CONSTS, c);
    
        backwarderror = (norm(K*c-f+h)/norm(-f+h));
        "Backward error: " + backwarderror
        forwardcu = norm(c(1:M)-c_exact(1:M))/norm(c_exact(1:M));
        forwardcv = norm(c(M:end)-c_exact(M:end))/norm(c_exact(M:end));
    
    
        it = it +1;
    end


    max_cu = max(abs(c(1:M)-c_exact(1:M))./c_exact(1:M));
    max_cv = max(abs(c(M+1:end)-c_exact(M+1:end))./c_exact(M+1:end));
    
    results(mesh_i, :) = [max_cu max_cv];

end

loglog(element_sizes, results(:, 1), '-or');
hold on;
loglog(element_sizes, results(:, 2), '-ob');
loglog(element_sizes, 1e-5*element_sizes.^2);

xlabel("Element size (edge length in mm)");
ylabel("Max relative error");
legend(["O2", "CO2", "Quadratic convergence"]);
