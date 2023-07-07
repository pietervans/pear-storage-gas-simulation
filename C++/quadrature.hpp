#ifndef quadrature
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <Eigen/Dense>
#define quadrature

/**
 * Evaluates the three basis functions N1, N2, N3 in the standard triangle
 * at each of the quadrature points of a certain quadrature rule.
 * The results is an N by 3 matrix, along with N weights, where the number of points
 * N is determined by the order.
 * % From "The Finite Element Method for Engineers" (2001, Wiley), p. 203
 */

inline std::pair<Eigen::MatrixXd, Eigen::ArrayXd> get_quadrature2d(int order)
{
    Eigen::MatrixX3d coords;
    Eigen::ArrayXd w;

    if (order == 1)
    {
        Eigen::Matrix<double, 1, 3> coords2;
        Eigen::Array<double, 1, 1> w2;
        coords2 << 1. / 3., 1. / 3., 1. / 3.;
        w2 << 1. / 2. * 1.;
        coords = coords2;
        w = w2;
    }
    else if (order == 2)
    {
        Eigen::Matrix<double, 3, 3> coords2;
        Eigen::Array<double, 3, 1> w2;
        coords2 << 1. / 2., 1. / 2., 0.,
            0, 1. / 2., 1. / 2.,
            1. / 2., 0., 1. / 2.;
        w2 << 1. / 3., 1. / 3., 1. / 3.;
        w2 *= 1. / 2.;
        coords = coords2;
        w = w2;
    }
    else if (order == 3)
    {
        Eigen::Matrix<double, 4, 3> coords2;
        Eigen::Array<double, 4, 1> w2;
        coords2 << 1. / 3., 1. / 3., 1. / 3.,
            0.6, 0.2, 0.2,
            0.2, 0.6, 0.2,
            0.2, 0.2, 0.6;
        w2 << -27. / 48., 25. / 48., 25. / 48., 25. / 48.;
        w2 *= 1. / 2.;
        coords = coords2;
        w = w2;
    }
    else
    {
        throw std::invalid_argument("Given order not supported!");
    }

    Eigen::Matrix<double, 3, 2> pts_standard;
    pts_standard << 0, 0, 1, 0, 0, 1;
    Eigen::MatrixX2d points = coords * pts_standard;

    // Evaluate points at Basis functions (cf. matlab)
    Eigen::MatrixXd N(w.rows(), 3);
    N(Eigen::all, 1) = points(Eigen::all, 0);
    N(Eigen::all, 2) = points(Eigen::all, 1);
    for (int i = 0; i < N.rows(); i++)
    {
        N(i, 0) = 1. - N(i, 1) - N(i, 2);
    }

    std::pair<Eigen::MatrixX3d, Eigen::ArrayXd> res(N, w);
    return res;
}

/**
 * Evaluates the two basis functions N1, N2 in the standard interval [0, 1]
 * at each of the quadrature points of a certain quadrature rule.
 * The results is an N by 2 matrix, along with N weights, where the number of points
 * N is determined by the order.
 * Uses Gaussian quadrature to get an integration rule that integrates all polynomials of degree
 * @param order or less exactly
 */
inline std::pair<Eigen::MatrixXd, Eigen::ArrayXd> get_quadrature1d(int order = 3)
{
    if (order >= 10)
    {
        throw std::invalid_argument("Given order not supported!");
    }

    int nb_of_points = (int)std::ceil((order + 1.) / 2.);
    // Gauss-Legendre quadrature
    Eigen::ArrayXd points(nb_of_points);
    Eigen::ArrayXd w(nb_of_points);
    double a, b, c, d, e;
    switch (nb_of_points)
    {
    case 1:
        points << 0.;
        w << 2.;
        break;
    case 2:
        a = 1. / std::sqrt(3.);
        points << -a, a;
        w << 1., 1.;
        break;
    case 3:
        a = std::sqrt(3. / 5.);
        points << -a, 0., a;
        w << 5. / 9., 8. / 9., 5. / 9.;
        break;
    case 4:
        a = 2. / 7. * std::sqrt(6. / 5.);
        b = std::sqrt(3. / 7. - a);
        c = std::sqrt(3. / 7. + a);
        d = (18 + std::sqrt(30.)) / 36.;
        e = (18 - std::sqrt(30.)) / 36.;
        points << -c, -b, b, c;
        w << e, d, d, e;
        break;
    case 5:
        a = 2. * std::sqrt(10. / 7.);
        b = 1. / 3. * std::sqrt(5. - a);
        c = 1. / 3. * std::sqrt(5. + a);
        d = (322. + 13. * std::sqrt(70.)) / 900.;
        e = (322. - 13. * std::sqrt(70.)) / 900.;
        points << -c, -b, 0., b, c;
        w << e, d, 128. / 225., d, e;
        break;
    }

    // Transform to [0, 1]
    points = points / 2. + 1. / 2.;
    w /= 2.;

    // Evaluate basis functions
    Eigen::MatrixXd N = Eigen::MatrixXd(nb_of_points, 2);
    N.col(0) = 1. - points;
    N.col(1) = points;

    std::pair<Eigen::MatrixX2d, Eigen::ArrayXd> res(N, w);
    return res;
}
#endif
