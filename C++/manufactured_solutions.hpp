#ifndef manufactured_solutions
#define manufactured_solutions

#include <functional>
#include <cmath>
#include "triangulation.hpp"
#include "constants.hpp"
#include "newton_iteration.hpp"
#include "quadrature.hpp"
#include "functions_ms.hpp"

template <typename F, typename F2>
class ManufacturedSolution {
    public:
        const F& Cu; //Chosen oxygen concentration
        const F& Cv; //Chosen Carbon dioxide concentration
        // Gradients
        const F& dCu_dr;
        const F& dCu_dz;
        const F& dCv_dr;
        const F& dCv_dz;
        const F& dCu_dr2;
        const F& dCu_dz2;
        const F& dCv_dr2;
        const F& dCv_dz2;
        const F2& Ru;
        const F2& Rv;
        const struct Constants& consts;

        inline ManufacturedSolution(
            F const& Cu_, F const& Cv_,
            F const& dCu_dr_, F const& dCu_dz_,
            F const& dCv_dr_, F const& dCv_dz_,
            F const& dCu_dr2_, F const& dCu_dz2_,
            F const& dCv_dr2_, F const& dCv_dz2_,
            F2 const& Ru_, F2 const& Rv_,
            struct Constants const& consts_
        ):
        Cu(Cu_), Cv(Cv_),
        dCu_dr(dCu_dr_), dCu_dz(dCu_dz_),
        dCv_dr(dCv_dr_), dCv_dz(dCv_dz_),
        dCu_dr2(dCu_dr2_), dCu_dz2(dCu_dz2_),
        dCv_dr2(dCv_dr2_), dCv_dz2(dCv_dz2_),
        Ru(Ru_), Rv(Rv_),
        consts(consts_)
        {};


        /**
         * Computes the f-term in the system of non-linear equations given the chosen manufactured solution.
         * There is a contribution from two terms:
         *      - The forcing term Q to make the PDE work.
         *      - The third integral over the boundary is now set by the chosen solution.
        */
        inline Eigen::VectorXd compute_forcing_term(
                triangulation const& tr, struct Constants const& consts, int order=3) {
            int M = tr.points.rows();
            int nb_of_triangles = tr.connectivity.rows();
            Eigen::VectorXd f(2*M);

            compute_prod_term_contribution(f, tr, consts, order);
            compute_boundary_contribution(f, tr, consts, order);
            
            return f;
        }

        inline Eigen::VectorXd get_c_exact(triangulation const& tr) {
            int M = tr.points.rows();
            Eigen::VectorXd c(2*M);
            for (int i=0; i<M; i++) {
                auto pt = tr.points(i,Eigen::all);
                c(i) =   Cu(pt(0),pt(1));
                c(M+i) = Cv(pt(0),pt(1));
            }
            return c;
        }

    private:
        //** Add contributions from integral over production term to forcing term f;
        inline void compute_prod_term_contribution(
                Eigen::VectorXd& f, triangulation const& tr, struct Constants const& consts, int order=3) {
            int M = tr.points.rows();
            int nb_of_triangles = tr.connectivity.rows();

            std::pair<Eigen::MatrixX3d, Eigen::ArrayXd> pair = get_quadrature2d(order);
            auto N = pair.first;
            auto w = pair.second;

            Eigen::ArrayXd qu(w.size());
            Eigen::ArrayXd qv(w.size());
        
            //Loop over triangles
            for (int k = 0; k < nb_of_triangles; k++) {
                Eigen::Array3i inds = tr.connectivity.row(k);
                Eigen::Matrix<double, 3, 2> rz = tr.points(inds, Eigen::all);
                Eigen::Vector3d r = rz.col(0);
                Eigen::Vector3d z = rz.col(1);
                double det_j = std::abs(det_jacobian(rz));

                Eigen::ArrayXd Nr = N*r; //Evaluate r at each quadrature point
                Eigen::ArrayXd Nz = N*z; //Evalute z at each quadrature point
                for (int i=0; i<w.size(); i++) { // Evaluate Q at each quadrature point
                    qu(i) = Qu(Nr(i), Nz(i));
                    qv(i) = Qv(Nr(i), Nz(i));
                }

                for (int m=0; m < 3; m++) {
                    int ind_m = inds(m);
                    Eigen::ArrayXd Nm = N.col(m);

                    f(ind_m) -= det_j*(qu*Nm*w).sum();
                    f(ind_m + M) -= det_j*(qv*Nm*w).sum();
                }
            }
        }

        /**
         * Add contributions from integral over boundary to forcing term f
        */
        inline void compute_boundary_contribution(
                Eigen::VectorXd& f, triangulation const& tr, struct Constants const& consts, int order=3) {

            int M = tr.points.rows();
            std::pair<Eigen::MatrixX2d, Eigen::ArrayXd> pair = get_quadrature1d(9);
            auto N = pair.first;
            auto w = pair.second;

            //Loop over boundary
            Eigen::MatrixX2i outer_boundary = tr.outer_boundary();
            for (int k=0; k< outer_boundary.rows(); k++) {
                Eigen::Array2i inds = outer_boundary.row(k);
                Eigen::Matrix<double, 2, 2> rz = tr.points(inds, Eigen::all);

                double deltar = rz(1, 0) - rz(0, 0);
                double deltaz = rz(1, 1) - rz(0, 1);
                double E2 = std::sqrt(std::pow(deltar,2) + std::pow(deltaz,2));

                //Compute outward normal vector
                double normalr = deltaz/E2;
                double normalz = -deltar/E2;

                Eigen::ArrayXd Nr = N*rz(Eigen::all, 0);
                Eigen::ArrayXd Nz = N*rz(Eigen::all, 1);
                Eigen::ArrayXd indepu(w.size()); //Terms independent from basis function in integrand
                Eigen::ArrayXd indepv(w.size()); //Terms independent from basis function in integrand
                for (int i=0; i<w.size(); i++) { //Evaluate term at each quadrature point
                    //Integral for oxygen
                    //evaluate gradient at quadrature point
                    double gradr = dCu_dr(Nr(i), Nz(i));
                    double gradz = dCu_dz(Nr(i), Nz(i));
                    indepu(i) = E2*Nr(i)*(normalr*consts.sigma_ur*gradr + normalz*consts.sigma_uz*gradz);
                    //Integral for carbon dioxide
                    //evaluate gradient at quadrature point
                    gradr = dCv_dr(Nr(i), Nz(i));
                    gradz = dCv_dz(Nr(i), Nz(i));
                    indepv(i) = E2*Nr(i)*(normalr*consts.sigma_vr*gradr + normalz*consts.sigma_vz*gradz);
                }

                for (int m=0; m<2; m++) {
                    int ind_m = inds(m);
                    Eigen::ArrayXd Nm = N.col(m);
                    f(ind_m) += (w*Nm*indepu).sum();
                    f(ind_m + M) += (w*Nm*indepv).sum();
                }
            }
        }

        /**
         * Production term in first PDE (right-hand side)
        */
        inline double Qu(double r, double z) {
            return consts.sigma_ur*(dCu_dr(r,z) + r*dCu_dr2(r,z))
               + consts.sigma_uz*(r*dCu_dz2(r,z))
               - r*Ru(r, z, consts);
        }

        /**
         * Production term in second PDE (right-hand side)
        */
        inline double Qv(double r, double z) {
            return consts.sigma_vr*(dCv_dr(r,z) + r*dCv_dr2(r,z))
               + consts.sigma_vz*(r*dCv_dz2(r,z))
               + r*Rv(r, z, consts);
        }

};

#endif