#include "triangulation.hpp"
#include <iostream>
#include "newton_iteration.hpp"
#include "fixed_terms.hpp"
#include "manufactured_solutions.hpp"
#include "functions_ms.hpp"
#include "io_operations.hpp"
#include <cmath>

/**
 * Return the adapted L2 "function norm", see section 5.6 from
 * https://l.messenger.com/l.php?u=http%3A%2F%2F160592857366.free.fr%2Fjoe%2Febooks%2FMechanical%2520Engineering%2520Books%2520Collection%2FFINITE%2520ELEMENT%2520ANALYSIS%2FA%2520first%2520corse%2520in%2520finite%2520element%2520analysis.pdf&h=AT3uAV82KJeaH8oH4LT_ESM6-tQ6MipKuOK9gnSD2Au-XSuSKODT5IgbAfHZSRQHJ9JfucxhdAs6r_hoqjdiTyQt6ZMEWMg-hW4bkJ9XN5X98c6Ei5mtd7LDmMSmSvw1SQ3X_dCi9tFpvg
*/
double norm_l2_adapted(Eigen::VectorXd v) {
    return std::sqrt(v.array().square().sum()/v.rows());
}

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
    
    ManufacturedSolution ms = ManufacturedSolution(
        Cu, Cv,
        dCu_dr,  dCu_dz,  dCv_dr,  dCv_dz,
        dCu_dr2, dCu_dz2, dCv_dr2, dCv_dz2,
        Ru, Rv, consts
    );

    std::list<std::string> mesh_list = std::list<std::string>();
    mesh_list.push_back("uniform_5mm");
    mesh_list.push_back("uniform_3mm");
    mesh_list.push_back("uniform_2mm");
    mesh_list.push_back("uniform_1mm");
    mesh_list.push_back("uniform_0p5mm");
    mesh_list.push_back("uniform_0p25mm");
    Eigen::ArrayXd element_sizes(mesh_list.size());
    element_sizes << 5.,3.,2.,1.,0.5,0.25;

    Eigen::MatrixXd results = Eigen::MatrixXd(mesh_list.size(), 3);
    results(Eigen::all, 0) = element_sizes;

    int i = 0;
    for (std::string mesh : mesh_list) {

        std::cout << "Finding solution for mesh " << mesh << std::endl;
        triangulation tr = get_triangulation(mesh);

        //Find solution
        int M = tr.points.rows();
        Eigen::SparseMatrix<double> J(2*M, 2*M);
        Eigen::SparseMatrix<double> K(2*M, 2*M);
        Eigen::VectorXd h = Eigen::VectorXd(2*M);

    
        integrals::fill_K_int1(K, tr, consts);
        Eigen::VectorXd f = ms.compute_forcing_term(tr, consts);
        Eigen::VectorXd c = integrals::get_initial_c_only_cu(K, f, tr, consts);
        newton_iter(c, K, f, tr, consts, 1.e-13, 10);

        
        //Take L-norm
        Eigen::VectorXd cu_exact(M);
        Eigen::VectorXd cv_exact(M);

        for (int j=0; j<M; j++) {
            Eigen::Array2d rz = tr.points(j, Eigen::all);
            double r = rz(0); double z = rz(1);

            cu_exact(j) = ms.Cu(r, z);
            cv_exact(j) = ms.Cv(r, z); 
        }
        
        results(i, 1) = norm_l2_adapted(c(Eigen::seq(0, M-1))-cu_exact)/norm_l2_adapted(cu_exact);
        results(i, 2) = norm_l2_adapted(c(Eigen::seq(M, 2*M-1))-cv_exact)/norm_l2_adapted(cv_exact);
        
        i++;
    }

    io::write_results(results, "../convergence_test.txt");
}
