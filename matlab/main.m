clear

CONSTS = get_constants('orchard');
tr = get_triangulation('adaptive_fine');
% tr = get_triangulation('uniform_5mm');
% tr = get_triangulation('uniform_3mm');
% tr = get_triangulation('uniform_2mm'); % NO CONVERGENCE!
% tr = get_triangulation('uniform_1mm'); % CONVERGENCE!
% tr = get_triangulation('uniform_0p5mm');
% tr = get_triangulation('uniform_0p25mm');

quadrature_order = 3;

tic()
M = size(tr.Points,1);
Ku = spalloc(M,M,7*M);
Kv = spalloc(M,M,7*M);
fu = zeros(M,1);
fv = zeros(M,1);

% Integral 1 in assignment (Equations (5) and (6))
Ku = Ku + get_K1(tr, CONSTS.sigma_ur, CONSTS.sigma_uz);
Kv = Kv + get_K1(tr, CONSTS.sigma_vr, CONSTS.sigma_vz);

% Integral 3 in assigment
[Ku_temp, fu_temp] = get_Kf3(tr, CONSTS.rho_u, CONSTS.C_uamb);
Ku = Ku + Ku_temp;
fu = fu + fu_temp;
[Kv_temp, fv_temp] = get_Kf3(tr, CONSTS.rho_v, CONSTS.C_vamb);
Kv = Kv + Kv_temp;
fv = fv + fv_temp;

c = getInitialC(Ku, Kv, fu, fv, tr, CONSTS);
"Get initial c: " + toc()

change = Inf;
while change > 1e-7
    
    % Linearization of nonlinear function
    tic();
    [J, h_ck] = get_Jacobian_2(quadrature_order, tr, CONSTS, c);
    "Get Jacobian: " + toc()
    
    K = [Ku, sparse(M,M); sparse(M,M) Kv];

    backwarderror = norm(K*c - [fu; fv] + h_ck)/norm(h_ck)

    % Solve system
    A = K + J;
    b = [fu; fv] - h_ck + J*c;
    
    %% Update c
    c_next = A\b;
    change = norm(c_next - c)/norm(c_next);
    c = c_next;

    
end

cu = c(1:M); cv = c(M+1:end);
subplot(1,2,1)
patch( ...
    'Faces',tr.ConnectivityList, ...
    'Vertices',tr.Points, ...
    'FaceVertexCData',cu, ...
    'FaceColor','interp' ...
)
colorbar
title('O2')
subplot(1,2,2)
patch( ...
    'Faces',tr.ConnectivityList, ...
    'Vertices',tr.Points, ...
    'FaceVertexCData',cv, ...
    'FaceColor','interp' ...
)
colorbar
title('CO2')

