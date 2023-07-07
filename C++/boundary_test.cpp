#include "triangulation.hpp"
#include <iostream>
#include "newton_iteration.hpp"
#include "fixed_terms.hpp"
#include "manufactured_solutions.hpp"
#include "io_operations.hpp"

int main(int argc, char *argv[]) {
    struct Constants consts;

    // Parse command-line arguments
    // Optional argument 1 specifies the constants
    if (argc == 1) {
        consts = Constants("orchard");
    } else if (argc == 2) {
        consts = Constants(argv[1]);
    } else {
        throw std::invalid_argument("Too many command-line arguments!");
    }

    std::list<std::string> mesh_list = std::list<std::string>();
    mesh_list.push_back("uniform_1mm");
    mesh_list.push_back("uniform_0p5mm");
    mesh_list.push_back("uniform_0p25mm");

    Eigen::MatrixXd results = Eigen::MatrixXd(mesh_list.size(), 3);


    int mesh_i = 0;
    for (std::string mesh : mesh_list) {

        triangulation tr = get_triangulation(mesh);
        int M = tr.points.rows();
        Eigen::SparseMatrix<double> J(2*M, 2*M);
        Eigen::SparseMatrix<double> K(2*M, 2*M);
        Eigen::VectorXd f = Eigen::VectorXd(M*2);
        integrals::fill_K_int1(K, tr, consts);
        integrals::fill_Kf_int3(K, f, tr, consts);
        Eigen::VectorXd c = integrals::get_initial_c(K, f, tr, consts);
        Eigen::VectorXd h = Eigen::VectorXd(M*2);
        newton_iter(c, K, f, tr, consts);

        
        Eigen::MatrixXi outer_boundary = tr.outer_boundary(true, true);
        Eigen::VectorXd cu_errors(outer_boundary.rows()), cv_errors(outer_boundary.rows());
        Eigen::VectorXd cu(3), cv(3);
        Eigen::Matrix<double, 3, 2> rz, dPhi;

        //Compute error on middle of boundary intervals
        for (int i = 0; i<outer_boundary.rows(); i++) {
            Eigen::Array3i inds = outer_boundary.row(i);
            rz = tr.points(inds, Eigen::all);

            double det_j = det_jacobian(rz);
            integrals::get_partial_derivatives(dPhi, rz, det_j);
            cu = c(inds);
            cv = c(inds+M);

            Eigen::Array2d derivatives_cu = dPhi.transpose()*cu;
            Eigen::Array2d derivatives_cv = dPhi.transpose()*cv;
            double deltar = rz(1, 0) - rz(0, 0);
            double deltaz = rz(1, 1) - rz(0, 1);
            double E2 = std::sqrt(std::pow(deltar,2) + std::pow(deltaz,2));


            //Compute outward normal vector
            double normalr = deltaz/E2;
            double normalz = -deltar/E2;

            //Right hand side and left hand side of boundary conditions
            double left_u = -(normalr*consts.sigma_ur*derivatives_cu(0) + normalz*consts.sigma_uz*derivatives_cu(1));
            double left_v = -(normalr*consts.sigma_vr*derivatives_cv(0) + normalz*consts.sigma_vz*derivatives_cv(1));

            double right_u = (deltar == 0.) ? 0. : consts.rho_u*((cu(0)+cu(1))/2. - consts.C_uamb);
            double right_v = (deltar == 0.) ? 0. : consts.rho_v*((cv(0)+cv(1))/2. - consts.C_vamb);

            cu_errors(i) = std::abs(left_u-right_u);
            cv_errors(i) = std::abs(left_v-right_v);


        };
        
        //Get maximum triangle size
        double max_size = 0;
        for (int k=0; k < tr.connectivity.rows(); k++) {
            Eigen::Array3i inds = tr.connectivity.row(k);
            Eigen::Matrix<double, 3, 2> rz = tr.points(inds, Eigen::all);
            double size = det_jacobian(rz)/2.;
            if (size > max_size) max_size = size;
        };
        results(mesh_i, 0) = max_size;

        results(mesh_i, 1) = cu_errors.maxCoeff();
        results(mesh_i, 2) = cv_errors.maxCoeff();

        mesh_i++;
    }
        
    io::write_results(results, "../boundary_test.txt");

    return 0;
}