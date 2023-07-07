#ifndef fixed_terms
#define fixed_terms

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include "triangulation.hpp"
#include <cassert>
#include "newton_iteration.hpp"

namespace integrals {

//Compute the partial derivatives of the basis functions with respect to r and z;
inline void get_partial_derivatives(Eigen::Matrix<double, 3, 2>& dPhi, Eigen::Matrix<double, 3, 2> const& rz, 
                                double det_j) {
    // Partial derivatives to r
    dPhi(0, 0) = (rz(1, 1)-rz(2, 1))/det_j;
    dPhi(1, 0) = (rz(2, 1)-rz(0, 1))/det_j;
    dPhi(2, 0) = -(rz(1, 1)-rz(0, 1))/det_j;

    //Parital derivatives to z
    dPhi(0, 1) = (rz(2, 0)-rz(1, 0))/det_j;
    dPhi(1, 1) = -(rz(2, 0)-rz(0, 0))/det_j;
    dPhi(2, 1) = (rz(1, 0)-rz(0, 0))/det_j;
}


inline double integral2d1(int i, int j, Eigen::Matrix<double, 3, 2> const& rz,
        Eigen::Matrix<double, 3, 2> const& dPhi, double det_j, double sigma_r, double sigma_z) {
    
    double partial_factor = sigma_r*dPhi(i, 0)*dPhi(j, 0) + sigma_z*dPhi(i,1)*dPhi(j,1);
    double res = rz(Eigen::all, 0).sum()/6.*partial_factor*det_j;
    return res;
}

inline double integral2d2(int i, int j, Eigen::Matrix<double, 3, 2> const& rz, double det_j) {
    if (i == j) {
        return 1./60.*(rz(Eigen::all,0).sum() + 2.*rz(i,0))*det_j;
    } else {
        return 1./120.*(rz(Eigen::all,0).sum() + rz(i,0) + rz(j,0))*det_j;
    }
}

inline double integral2d3(int j, Eigen::Matrix<double, 3, 2> const& rz, double det_j) {
    return (rz(Eigen::all,0).sum() + rz(j,0))*det_j/24.;
}

inline double integral1d1(int j, Eigen::Matrix<double, 2, 2> const& rz, double norm_E) {
    return (rz(0,0) + rz(1,0) + rz(j,0))/6.*norm_E;
}

inline double integral1d2(int i, int j, Eigen::Matrix<double, 2, 2> const& rz, double norm_E) {
    Eigen::Vector2d r = rz.col(0);
    if (i == j) {
        return 1./12.*(r(0)+r(1)+2*r(i))*norm_E;
    } else {
        return 1./12.*(r(0)+r(1))*norm_E;
    }
}


inline void fill_K_int1(
        Eigen::SparseMatrix<double>& K, triangulation const& Tr, 
        struct Constants const& consts) {
    int M = Tr.points.rows();
    int nb_of_triangles = Tr.connectivity.rows();
    assert(K.rows() == M*2 && K.cols() == M*2);
    Eigen::SparseMatrix<double> deltaK(2*M, 2*M);
    std::vector<Triplet> tripletList;

    for (int k=0; k < nb_of_triangles; k++) {
        Eigen::Array3i inds = Tr.connectivity.row(k);
        Eigen::Matrix<double, 3, 2> rz = Tr.points(inds, Eigen::all);
        Eigen::Matrix<double, 3, 2> dPhi;

        double det_j = det_jacobian(rz);
        get_partial_derivatives(dPhi, rz, det_j);

        for (int i=0; i<3; i++) {
            int ind_i = inds(i);
            for (int j=0; j<3; j++) {
                int ind_j = inds(j);
                double resU = integral2d1(i, j, rz, dPhi, det_j, consts.sigma_ur, consts.sigma_uz);
                double resV = integral2d1(i, j, rz, dPhi, det_j, consts.sigma_vr, consts.sigma_vz);
                
                tripletList.push_back(Triplet(ind_j, ind_i, resU));
                tripletList.push_back(Triplet(ind_j + M, ind_i + M, resV));
            }
        }
    }
    deltaK.setFromTriplets(tripletList.begin(), tripletList.end());
    K += deltaK;
}

/**
 * K and f passed by reference! Make sure not to change them for Newton iterations!
*/
inline void fill_Kf_int2(
        Eigen::SparseMatrix<double>& K, Eigen::VectorXd& f,
        triangulation const& Tr, struct Constants const& consts) {
    int M = Tr.points.rows();
    assert(K.rows() == M*2 && K.cols() == M*2);
    assert(f.rows() == M*2);
    int nb_of_triangles = Tr.connectivity.rows();

    Eigen::SparseMatrix<double> deltaK(2*M, 2*M);
    std::vector<Triplet> tripletList;

    for (int k=0; k < nb_of_triangles; k++) {
        Eigen::Array3i inds = Tr.connectivity.row(k);
        Eigen::Matrix<double, 3, 2> rz = Tr.points(inds, Eigen::all);
        double det_j = det_jacobian(rz);
        for (int j=0; j<3; j++) {
            int ind_j = inds(j);
            double integral_2d3 = integral2d3(j,rz,det_j);
            f(M+ind_j) += consts.V_mfv*integral_2d3;
            for (int i=0; i<3; i++) {
                int ind_i = inds(i);
                double integral_2d2 = integral2d2(i,j,rz,det_j);

                tripletList.push_back(Triplet(ind_j,ind_i, consts.V_mu/consts.K_mu*integral_2d2));
                // Update lower left block of K:
                tripletList.push_back(Triplet(M+ind_j,ind_i, - consts.r_q*consts.V_mu/consts.K_mu*integral_2d2));

            }
        }
    }
    deltaK.setFromTriplets(tripletList.begin(), tripletList.end());
    K += deltaK;
}
    

inline void fill_Kf_int3(
        Eigen::SparseMatrix<double>& K, Eigen::VectorXd& f,
        triangulation const& Tr, struct Constants const& consts) {
    int M = Tr.points.rows();
    assert(K.rows() == M*2 && K.cols() == M*2);
    assert(f.rows() == M*2);
    Eigen::MatrixX2i boundary_inds = Tr.outer_boundary();
    int nb_of_edges = boundary_inds.rows();
    Eigen::SparseMatrix<double> deltaK(2*M, 2*M);
    std::vector<Triplet> tripletList;


    for (int k=0; k < nb_of_edges; k++) {
        Eigen::Array2i inds = boundary_inds.row(k);
        Eigen::Matrix<double, 2, 2> rz = Tr.points(inds, Eigen::all);
        double norm_E = (rz(0, Eigen::all) - rz(1,Eigen::all)).norm();
        for (int j=0; j < 2; j++) {
            int ind_j = inds(j);
            double integral_j = integral1d1(j,rz,norm_E);
            f(ind_j)   += consts.rho_u*consts.C_uamb*integral_j;
            f(M+ind_j) += consts.rho_v*consts.C_vamb*integral_j;
            for (int i=0; i < 2; i++) {
                int ind_i = inds(i);
                double integral_ij = integral1d2(i,j,rz,norm_E);

                tripletList.push_back(Triplet(ind_j,ind_i, consts.rho_u*integral_ij));
                tripletList.push_back(Triplet(M+ind_j,M+ind_i,consts.rho_v*integral_ij));
            }
        }
    }
    deltaK.setFromTriplets(tripletList.begin(), tripletList.end());
    K += deltaK;
}

/**
 * K and f are passed by value because they should not be changed for caller!
 * K and f should contain contributions from integral 1 and 3 (cf. matlab)
*/
inline Eigen::VectorXd get_initial_c(Eigen::SparseMatrix<double> K, Eigen::VectorXd f,
        triangulation const& Tr, struct Constants const& consts) {
    int M = Tr.points.rows();
    assert(K.rows() == M*2 && K.cols() == M*2);
    assert(f.rows() == M*2);
    fill_Kf_int2(K, f, Tr, consts);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(K);
    solver.factorize(K);
    return solver.solve(f);
}

/**
 * Return initial vector c where only cu is computed, cv is set to zero.
 * This avoids issues caused by ill-conditioned matrix Kv.
*/
inline Eigen::VectorXd get_initial_c_only_cu(Eigen::SparseMatrix<double> K, Eigen::VectorXd f,
        triangulation const& Tr, struct Constants const& consts) {
    int M = Tr.points.rows();
    fill_Kf_int2(K, f, Tr, consts);

    // Extract upper left block of K into K1
    std::vector<Triplet> tripletList;
    tripletList.reserve(7*M);
    for (int k=0; k<K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K,k); it; ++it) {
            int ind_row = it.row();
            int ind_col = it.col();
            if (ind_row < M && ind_col < M) {
                tripletList.push_back(Triplet(ind_row, ind_col, it.value()));
            }
        }
    }
    Eigen::SparseMatrix<double> K1(M,M);
    K1.setFromTriplets(tripletList.begin(), tripletList.end());

    // Solve for Cu
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(K1);
    solver.factorize(K1);
    auto fu = f(Eigen::seq(0,M-1));
    auto cu = solver.solve(fu);
    Eigen::VectorXd c(2*M);
    c.setZero();
    c(Eigen::seq(0,M-1)) = cu;
    return c;
}

}

#endif
