#ifndef io_operations_hpp
#define io_operations_hpp
#include <Eigen/Dense>
#define MAXBUFSIZE  ((int) 1e8)
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
namespace io {

/**
 * Read matrix of arbitrary dimension from .txt file and return corresponding Eigen::Matrix object.
 * Adapted from https://gist.github.com/zodsoft/52c738e60c264bdb2a5967b5ca0d8de6
*/
inline Eigen::MatrixXd read_matrix(std::string filename) {
    std::ifstream in(filename);
    std::string line;
    int row = 0;
    int col = 0;
    int nb_cols = 0;
    double* buffer = new double[MAXBUFSIZE];

    while (std::getline(in, line)) {
        char *ptr = (char *) line.c_str();
        col = 0;
        char *start = ptr;
        for (int i = 0; i < line.length(); i++) {
            if (ptr[i] == ',') {
                if (nb_cols*row + col >= MAXBUFSIZE) {
                    throw std::runtime_error("Error: MAXBUFSIZE too small!");
                }
                buffer[ nb_cols*row + col++ ] = atof(start);
                start = ptr + i + 1;
            }
        }
        if (col==0) continue;
        if (nb_cols==0) nb_cols=col+1;
        buffer[ nb_cols*row + col ] = atof(start);
        row++;
    }
    in.close();
    int nb_rows = row;

    Eigen::MatrixXd res(nb_rows, nb_cols);
    for (int i = 0; i < nb_rows; i++) {
        for (int j = 0; j < nb_cols; j++) {
            res(i,j) = buffer[ nb_cols*i+j ];
        }
    }
    delete[] buffer;
    return res;
}

inline void write_results(Eigen::VectorXd const& c, std::string filename="../c_solution.txt") {
    std::ofstream myfile;
    myfile.open(filename);
    for (int i=0; i<c.rows(); i++) {
        myfile << c(i) << "\n";
    }
    myfile.close();
}

inline void write_results(Eigen::MatrixXd const& A, std::string filename="../convergence_test.txt") {
    std::ofstream myfile;
    myfile.open(filename);
    for (int i=0; i<A.rows(); i++) {
        for (int j=0; j<A.cols()-1; j++) {
            myfile << A(i,j) << ",";
        }
        myfile << A(i, A.cols()-1) << "\n";
    }
    myfile.close();
}


}

#endif