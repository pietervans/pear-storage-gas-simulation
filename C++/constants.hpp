#include <string>
#include <stdexcept>
#include <cmath>
#ifndef constants
#define constants


    struct Constants
    {
        double T_cel;
        double eta_u;
        double eta_v;
        double T;

        double sigma_ur = 2.8*1e-10;
        double sigma_uz = 1.1*1e-9;
        double sigma_vr = 2.32*1e-9;
        double sigma_vz = 6.97*1e-9;

        double R_g = 8.314;
        double T_ref = 293.15;
        double V_muref = 2.39*1e-4;
        double E_avmuref = 80200;
        double V_mu;

        double V_mfvref = 1.61*1e-4;
        double E_avmfvref = 56700.;
        double V_mfv;

        double K_mu = 0.4103;
        double K_mv = 27.2438;
        double K_mfu = 0.1149;
        double r_q = 0.97;
        double rho_u = 7*1e-7;
        double rho_v = 7.5*1e-7;
        double p_atm = 101300.;

        double C_uamb;
        double C_vamb;


        inline Constants(std::string storage_description)
        { 
            if ( storage_description == "orchard") {
                T_cel = 25.;
                eta_u = 20.8/100.;
                eta_v = 0.04/100.;
            } else if ( storage_description == "shelf_life") {
                T_cel = 20.;
                eta_u = 20.8/100.;
                eta_v = 0./100.;
            } else if ( storage_description == "refrigerator"){
                T_cel = 7.;
                eta_u = 20.8/100.;
                eta_v = 0./100.;
            } else if ( storage_description == "precooling") {
                T_cel = -1.;
                eta_u = 20.8/100.;
                eta_v = 0./100.;
            } else if ( storage_description == "disorder_inducing") {
                T_cel = -1.;
                eta_u = 2./100.;
                eta_v = 5./100.;
            } else if ( storage_description == "optimal_ca") {
                T_cel = -1.;
                eta_u = 2./100.;
                eta_v = 0.7/100.;
            } else {
                throw std::invalid_argument("Invalid storage description!");
            };

            T = T_cel + 273.15;
            V_mu = V_muref*std::exp(E_avmuref/R_g*(1./T_ref - 1./T));
            V_mfv = V_mfvref*exp(E_avmfvref/R_g*(1./T_ref - 1./T));

            C_uamb = p_atm*eta_u/(R_g*T);
            C_vamb = p_atm*eta_v/(R_g*T);
        };

        inline Constants() {}

        Constants& operator=(const Constants& c)
        {
            sigma_ur = c.sigma_ur;
            sigma_uz = c.sigma_uz;
            sigma_vr = c.sigma_vr;
            sigma_vz = c.sigma_vz;
            K_mu = c.K_mu;
            V_mu = c.V_mu;
            V_mfv = c.V_mfv;
            K_mv = c.K_mv;
            K_mfu = c.K_mfu;
            r_q = c.r_q;
            rho_u = c.rho_u;
            rho_v = c.rho_v;
            C_uamb = c.C_uamb;
            C_vamb = c.C_vamb;
            return *this;
        }
    };

#endif