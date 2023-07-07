#ifndef funs_ms
#define funs_ms

#include <cmath>
#include "constants.hpp"
using namespace std;

inline double Cu(double r, double z) {
    return 0.05 + 0.05*sin(10*r) + pow(z,2);
}
inline double Cv(double r, double z) {
    return 2 + 3*pow(r,2) + sin(15*z);;
}
inline double dCu_dr(double r, double z) {
    return 0.05*10*cos(10*r);
}
inline double dCu_dz(double r, double z) {
    return 2*z;
}
inline double dCv_dr(double r, double z) {
    return 6*r;
}
inline double dCv_dz(double r, double z) {
    return 15*cos(15*z);
}
inline double dCu_dr2(double r, double z) {
    return -10*10*0.05*sin(10*r);
}
inline double dCu_dz2(double r, double z) {
    return 2;
}
inline double dCv_dr2(double r, double z) {
    return 6;
}
inline double dCv_dz2(double r, double z) {
    return -15*15*sin(15*z);
}

inline double Ru(double r, double z, struct Constants const& consts) {
    return consts.V_mu*Cu(r, z)/(consts.K_mu + Cu(r,z))/(1. + Cv(r,z)/consts.K_mv);
}
inline double Rv(double r, double z, struct Constants const& consts) {
    return consts.r_q*Ru(r,z,consts) + consts.V_mfv/(1. + Cu(r,z)/consts.K_mfu);
}

#endif
