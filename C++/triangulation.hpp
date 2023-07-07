#ifndef triangulation_
#define triangulation_

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
#include <Eigen/Dense>
#include "io_operations.hpp"
#include <stdexcept>
#include <algorithm>
#include <unordered_set>
#include <list>
#define MAXBUFSIZE  ((int) 1e8)

/**
 * Object used for 2D triangle mesh
*/
class triangulation {
    public:
    inline triangulation(
        Eigen::MatrixX2d points_,
        Eigen::MatrixX3d connectivity_
    )
    : points(points_), connectivity(connectivity_.cast<int>())
    {
        //Make sure triangles are specified with points in counter-clockwise order
        for (int k=0; k < connectivity.rows(); k++) {
            Eigen::Array3i inds = connectivity.row(k);
            Eigen::Matrix<double, 3, 2> rz = points(inds, Eigen::all);
            double dr21 = rz(1, 0) - rz(0, 0);
            double dz21 = rz(1, 1) - rz(0, 1);
            double dr31 = rz(2, 0) - rz(0, 0);
            double dz31 = rz(2, 1) - rz(0, 1);
            
            double normal_direction = dr21*dz31 - dr31*dz21;
            if (normal_direction < 0) { // Flip triangle
                connectivity(k, 1) = inds(2);
                connectivity(k, 2) = inds(1);
            }
            
        }
    }

    inline triangulation() {}

    inline triangulation& operator=(const triangulation& t)
    {
        points = t.points;
        connectivity = t.connectivity;
        return *this;
    }

    Eigen::MatrixX2d points;
    Eigen::MatrixX3i connectivity;


    struct edgeHashFunction {
        size_t operator()(const std::pair<int,int>& edge) const {
            // Not perfectly collision resistant, but hopefully good enough
            return edge.first^edge.second; // bitwise XOR
        }
    };

    inline Eigen::MatrixXi outer_boundary(bool inlcude_gamma1 = false, bool add_third_index=false) const {
        // Return a matrix with rows containing the point-ids on the outer boundary of the pear
        std::unordered_map<std::pair<int,int>, std::tuple<int, int, int>,edgeHashFunction> edges;
        for (int k=0; k<connectivity.rows(); k++) {
            for (int i=0; i<2; i++) {
                int id_i = connectivity(k,i);
                for (int j=i+1; j<3; j++) {
                    int id_j = connectivity(k,j);
                    int id_third = connectivity(k, 3-i-j);
                    // Always store edges with ids in same order (counter-clockwise rotation)
                    std::pair edge_key(std::min(id_i,id_j),std::max(id_i,id_j)); //Key
                    
                    std::tuple edge(id_i, id_j, id_third);
                    if (i==0 && j==2) {
                        edge = std::tuple(id_j, id_i, id_third);
                    }//Flip.
                    

                    auto ind = edges.find(edge_key);
                    bool found = (ind != edges.end());
                    if (found) {
                        // Not on boundary! (2 triangles share this edge)
                        edges.erase(ind);
                    } else {
                        edges.insert(std::pair(edge_key, edge));
                    }
                }
            }
        }
        std::list<std::tuple<int,int, int>> outer_edges;
        for (auto const& edge : edges) {
            double r1 = points(std::get<0>(edge.second),0);
            double r2 = points(std::get<1>(edge.second),0);
            if (inlcude_gamma1 || (r1!=0 || r2!=0)) {
                outer_edges.push_back(edge.second);
            }
        }
        Eigen::MatrixXi boundary_edges(outer_edges.size(),add_third_index ? 3 : 2);
        int i = 0;
        for (auto const& edge : outer_edges) {
            boundary_edges(i,0) = std::get<0>(edge);
            boundary_edges(i,1) = std::get<1>(edge);
            if (add_third_index) {
                boundary_edges(i,2) = std::get<2>(edge);
            };
            i++;
        }
        
            
        return boundary_edges;
    }

};

/**
 * Read the triangulation object of the given type.
 * type == 'adaptive_[fine|rough]' or type == 'uniform_[5|3|1|0p5|0p25]mm']
*/
inline triangulation get_triangulation(std::string type) {
    std::string prefix = "../matlab/pear_meshes/data/";
    auto points =       io::read_matrix(prefix + "points_" +    type + ".txt");
    auto connectivity = io::read_matrix(prefix + "triangles_" + type + ".txt");
    // Convert to array since Eigen does not support scalar subtraction on Matrix
    Eigen::ArrayX3d connectivity_arr = connectivity;
    // Make sure indices start at 0 (as opposed to matlab)
    connectivity_arr -= 1;
    connectivity = connectivity_arr;
    triangulation tr(points,connectivity);
    return tr;
}

#endif
