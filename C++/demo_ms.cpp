#include "triangulation.hpp"
#include <iostream>
#include "newton_iteration.hpp"
#include "fixed_terms.hpp"
#include "manufactured_solutions.hpp"
#include "functions_ms.hpp"
#include "io_operations.hpp"

int main(int argc, char *argv[]) {
    struct Constants consts;
    triangulation tr;

    // Parse command-line arguments
    // Optional argument 1 is type of triangulation
    // Optional argument 2 specifies the constants
    if (argc == 1) {
        tr = get_triangulation("uniform_1mm");
        consts = Constants("orchard");
    } else if (argc == 2) {
        tr = get_triangulation(argv[1]);
        consts = Constants("orchard");
    } else if (argc == 3) {
        tr = get_triangulation(argv[1]);
        consts = Constants(argv[2]);
    } else {
        throw std::invalid_argument("Too many command-line arguments!");
    }
    
    int M = tr.points.rows();
    Eigen::SparseMatrix<double> J(2*M, 2*M);
    Eigen::SparseMatrix<double> K(2*M, 2*M);
    Eigen::VectorXd h = Eigen::VectorXd(2*M);

    ManufacturedSolution ms = ManufacturedSolution(
        Cu, Cv,
        dCu_dr,  dCu_dz,  dCv_dr,  dCv_dz,
        dCu_dr2, dCu_dz2, dCv_dr2, dCv_dz2,
        Ru, Rv, consts
    );

    std::cout << "Int 1" << std::endl;
    integrals::fill_K_int1(K, tr, consts);
    std::cout << "Compute forcing term" << std::endl;
    Eigen::VectorXd f = ms.compute_forcing_term(tr, consts);
    std::cout << "Compute initial c" << std::endl;
    Eigen::VectorXd c = integrals::get_initial_c_only_cu(K, f, tr, consts);
    std::cout << "Newton iterations" << std::endl;
    newton_iter(c, K, f, tr, consts);

    Eigen::VectorXd c_exact = ms.get_c_exact(tr);
    std::cout << "Average residue equals " << (c-c_exact).cwiseAbs().sum()/(2*M) << std::endl;

    io::write_results(c, "../c_solution.txt");
}
