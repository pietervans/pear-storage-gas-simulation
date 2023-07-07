#ifndef newton_iteration
#define newton_iteration

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include "triangulation.hpp"
#include "constants.hpp"
#include <Eigen/Dense>
#include <string>
#include <Eigen/SparseCore>
#include <cassert>
#include <iostream>
#include <Eigen/SparseLU>
#include "quadrature.hpp"
#include <chrono>
#include <fstream>


inline double det_jacobian(Eigen::Matrix<double, 3, 2> rz) {
        double diag = (rz(1, 0)-rz(0, 0))*(rz(2, 1)-rz(0, 1)) ;
        double diag2 = (rz(1, 1) - rz(0, 1))*(rz(2, 0) - rz(0, 0));
        return std::abs(diag - diag2);
}


typedef Eigen::Triplet<double> Triplet;
inline void compute_jacobian(Eigen::SparseMatrix<double>& J,Eigen::VectorXd const& c, 
        Eigen::VectorXd& h, triangulation const& Tr, struct Constants const& consts, int ord=3) {
    
    int M = Tr.points.rows();
    int nb_of_triangles = Tr.connectivity.rows();
    assert(c.rows() == M*2 && h.rows() == M*2);
    std::vector<Triplet> tripletList;
    tripletList.reserve(4*7*M);

    h.setZero();

    std::pair<Eigen::MatrixX3d, Eigen::ArrayXd> pair = get_quadrature2d(ord);
    auto N = pair.first;
    auto w = pair.second;

    Eigen::ArrayXd ones = Eigen::ArrayXd::Ones(w.size());

    for (int k=0; k < nb_of_triangles; k++) {
        Eigen::Array3i inds = Tr.connectivity.row(k);
        Eigen::Matrix<double, 3, 2> rz = Tr.points(inds, Eigen::all);
        Eigen::Vector3d r = rz.col(0);
        double det_j = det_jacobian(rz);


        Eigen::Vector3d cu = c(inds);
        Eigen::Vector3d cv = c(inds + M);

        Eigen::ArrayXd Nr =  N*r;
        Eigen::ArrayXd Ncu = N*cu;
        Eigen::ArrayXd Ncv = N*cv;

        Eigen::ArrayXd noemer1 = 1. + (1./consts.K_mv)*Ncv;
        Eigen::ArrayXd noemer2 = consts.K_mu + Ncu;
        Eigen::ArrayXd noemer3 = (1. + (1./consts.K_mfu)*Ncu).square();


        for (int m = 0; m<3; m++) {
            int ind_m = inds(m);
            Eigen::ArrayXd Nm = N.col(m);

            auto delta_integrand = Nr*Ncu*Nm/noemer1/noemer2;
            double delta = consts.V_mu*det_j*(w*delta_integrand).sum();

            auto epsilon_integrand = Nr*Nm/(1.+1./consts.K_mfu*Ncu);
            double epsilon = consts.V_mfv*det_j*(w*epsilon_integrand).sum();

            h(ind_m) += delta;
            h(ind_m + M) += -consts.r_q*delta - epsilon;

            for (int i=m; i<3; i++) {
                double ind_i = inds(i);
                Eigen::ArrayXd Ni = N.col(i);

                auto tripleN = Nr*Ni*Nm;
                
                auto alpha_integrand = tripleN/noemer1/(noemer2.square());
                double alpha = consts.V_mu*consts.K_mu*det_j*(w*alpha_integrand).sum();

                auto beta_integrand = Ncu*tripleN/(noemer1.square())/noemer2;
                double beta = -consts.V_mu/consts.K_mv*det_j*(w*beta_integrand).sum();

                auto gamma_integrand = tripleN/noemer3;
                double gamma = consts.V_mfv/consts.K_mfu*det_j*(w*gamma_integrand).sum();
                
                tripletList.push_back(Triplet(ind_m, ind_i, alpha));
                tripletList.push_back(Triplet(ind_m, ind_i+M, beta));
                tripletList.push_back(Triplet(ind_m+M, ind_i, -consts.r_q*alpha + gamma));
                tripletList.push_back(Triplet(ind_m+M, ind_i+M, -consts.r_q*beta));

                if (m != i) {
                    tripletList.push_back(Triplet(ind_i, ind_m, alpha));
                    tripletList.push_back(Triplet(ind_i, ind_m+M, beta));
                    tripletList.push_back(Triplet(ind_i+M, ind_m, -consts.r_q*alpha + gamma));
                    tripletList.push_back(Triplet(ind_i+M, ind_m+M, -consts.r_q*beta));
                }
            }
        }    
    }

    J.setFromTriplets(tripletList.begin(), tripletList.end());
    
}

inline void newton_iter(Eigen::VectorXd& c, Eigen::SparseMatrix<double> const& K, Eigen::VectorXd& f, 
                            triangulation const& Tr, struct Constants const& consts,
                            double treshhold=1e-11, int max_iterations=20) {
    int M = Tr.points.rows();
    assert(f.rows() == 2*M);
    assert(K.rows() == 2*M && K.cols() == 2*M);

    
    Eigen::VectorXd h(2*M);
    Eigen::VectorXd b(2*M);
    Eigen::SparseMatrix<double> A(2*M, 2*M);
    Eigen::SparseMatrix<double> J(2*M, 2*M);

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    int i = 0;
    double error = 100;
    while (error > treshhold && i < max_iterations) {

        //Compute jacobian
        compute_jacobian(J, c, h, Tr, consts);
        
        error = (K*c - f + h).norm()/(-f+h).norm();
        std::cout << "Error: " << error << std::endl;

        //Solve system
        b = f - h + J*c;
        A = K+J;
        
        if (i==0) solver.analyzePattern(A); //Pattern stays the same
        solver.factorize(A);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Was unable to factorize sparse matrix to solve system.");
        }
        
        c = solver.solve(b);
        i++;
    }
}


#endif
