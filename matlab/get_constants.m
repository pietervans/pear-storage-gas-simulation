function CONSTS = get_constants(storage_description)
switch storage_description
    case "orchard"
        T_cel = 25;
        eta_u = 20.8/100;
        eta_v = 0.04/100;
    case "shelf_life"
        T_cel = 20;
        eta_u = 20.8/100;
        eta_v = 0/100;
    case "refrigerator"
        T_cel = 7;
        eta_u = 20.8/100;
        eta_v = 0/100;
    case "precooling"
        T_cel = -1;
        eta_u = 20.8/100;
        eta_v = 0/100;
    case "disorder_inducing"
        T_cel = -1;
        eta_u = 2/100;
        eta_v = 5/100;
    case "optimal_ca"
        T_cel = -1;
        eta_u = 2/100;
        eta_v = 0.7/100;
    otherwise
        error("Invalid storage description!")
end
T = T_cel + 273.15;

sigma_ur = 2.8*1e-10;
sigma_uz = 1.1*1e-9;
sigma_vr = 2.32*1e-9;
sigma_vz = 6.97*1e-9;

R_g = 8.314;
T_ref = 293.15;
V_muref = 2.39*1e-4;
E_avmuref = 80200;
V_mu = V_muref*exp(E_avmuref/R_g*(1/T_ref - 1/T));

V_mfvref = 1.61*1e-4;
E_avmfvref = 56700;
V_mfv = V_mfvref*exp(E_avmfvref/R_g*(1/T_ref - 1/T));

K_mu = 0.4103;
K_mv = 27.2438;
K_mfu = 0.1149;
r_q = 0.97;
rho_u = 7*1e-7;
rho_v = 7.5*1e-7;
p_atm = 101300;

C_uamb = p_atm*eta_u/(R_g*T);
C_vamb = p_atm*eta_v/(R_g*T);

CONSTS = struct();
CONSTS.sigma_ur = sigma_ur;
CONSTS.sigma_uz = sigma_uz;
CONSTS.sigma_vr = sigma_vr;
CONSTS.sigma_vz = sigma_vz;
CONSTS.K_mu = K_mu;
CONSTS.V_mu = V_mu;
CONSTS.V_mfv = V_mfv;
CONSTS.K_mv = K_mv;
CONSTS.K_mfu = K_mfu;
CONSTS.r_q = r_q;
CONSTS.rho_u = rho_u;
CONSTS.rho_v = rho_v;
CONSTS.C_uamb = C_uamb;
CONSTS.C_vamb = C_vamb;
end
