#include "constants.hpp"
#include "triangulation.hpp"
#include <assert.h> 
#include <random>
#include "quadrature.hpp"
#include <iostream>

void test_constants() {
    Constants consts = Constants("optimal_ca");
    assert(consts.T == 272.15);
}

void test_triangulation() {
    triangulation tr = get_triangulation("uniform_5mm");
    auto bound = tr.outer_boundary(false); //No gamma 1
    assert(bound.rows() == 71-26);
    bound = tr.outer_boundary(true);
    assert(bound.rows() == 71);

    tr = get_triangulation("uniform_2mm");
    bound = tr.outer_boundary(false); //No gamma 1
    assert(bound.rows() == 140-60);
    bound = tr.outer_boundary(true);
    assert(bound.rows() == 140);
}

void test_1d_quadrature() {
    for (int degree = 0; degree < 10; degree++) {
        std::pair<Eigen::MatrixX2d, Eigen::ArrayXd> pair = get_quadrature1d(degree);
        auto N = pair.first;
        auto w = pair.second;
        
        //Evaluate random polynomial with points in second basis function f(xi) = xi;
        std::uniform_real_distribution<double> unif(0, 1);
        std::default_random_engine re;
        double int_exact = 0.;
        double int_quadrature = 0.;
        
        for (int i = 0; i < degree+1; i++) {
            double coeff = unif(re);
            int_exact += coeff/((double)i+1.);
            int_quadrature += (w*coeff*(N.col(1).array().pow(i))).sum();
        }
        
        assert(std::abs(int_exact - int_quadrature)/std::abs(int_exact) <= std::numeric_limits<double>::epsilon());
        

        //Make sure column 1 is first basis function f(xi) = 1-xi
        Eigen::ArrayXd difference = ((1.0 - N.col(1).array()) - N.col(0).array()).abs();
        Eigen::ArrayXd rel_difference = difference/N.col(1).array().abs();
        
        assert((rel_difference <= std::numeric_limits<double>::epsilon()).all());
    }
}

int main() {
    test_constants();
    test_triangulation();
    test_1d_quadrature();
    return 0;
}
