#include "triangulation.hpp"
#include <iostream>
#include "newton_iteration.hpp"
#include "fixed_terms.hpp"
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
    
    std::cout << "Computing initial concentration..." << std::endl;
    int M = tr.points.rows();
    Eigen::SparseMatrix<double> J(2*M, 2*M);
    Eigen::SparseMatrix<double> K(2*M, 2*M);
    Eigen::VectorXd f = Eigen::VectorXd(M*2);
    integrals::fill_K_int1(K, tr, consts);
    
    integrals::fill_Kf_int3(K, f, tr, consts);
    
    
    Eigen::VectorXd c = integrals::get_initial_c(K, f, tr, consts);
    Eigen::VectorXd h = Eigen::VectorXd(M*2);
    
    std::cout << "Starting newton iteration..." << std::endl;
    newton_iter(c, K, f, tr, consts);

    std::cout << "Saving results to c_solution.txt..." << std::endl;
    io::write_results(c, "../c_solution.txt");

    return 0;
}
