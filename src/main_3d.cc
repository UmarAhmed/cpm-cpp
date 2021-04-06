#include <armadillo>
#include <vector>
#include <iostream>
#include <algorithm>

#include "point_cloud.h"
#include "cpm_util.h"


using namespace arma;


int main() {
    // point cloud sphere
    constexpr int M = 50;
    constexpr double pi = 3.14159;
    
    const vec theta = linspace(0, 2 * pi, 2 * M);
    const vec phi = linspace(0, pi, M);

    std::vector<vec> surface;

    for (const auto t: theta) {
        for (const auto p: phi) {
            const vec q {sin(p) * sin(t), sin(p) * cos(t), cos(p) }; 
            surface.emplace_back(q);
        }
    }
    std::cout << "Number of points in surface: " << surface.size() << std::endl; 

    // create grid
    constexpr double dx = 0.1;
    constexpr double bound = 2;
    const vec pts = regspace(-bound, dx, bound);
    const int n = pts.size();
    const int N = n * n * n;
    
    std::cout << "Number of points in grid: " << N << std::endl; 

    // Compute warm start (list of surface pts within the bandwidth for each grid pt)
    const double threshold = 0.41231469; // TODO hard coded for now
    const std::vector<std::vector<int>> starters = warm_start_3d(pts, surface, threshold);


    // Check average starter size
    int mm = 0;
    for (int i = 0; i < starters.size(); i++) {
        mm += starters[i].size();
    }
    std::cout << "avg starter size: " << mm / N << std::endl;

    std::vector<vec> cp_pts;
    std::vector<int> band;

    for (int idx = 0; idx < N; idx++) {
        if (starters[idx].size() == 0) {
            continue;
        }

        const int i = idx % n;
        const int j = (idx / n) % n;
        const int k = idx / (n * n);

        const vec p {pts(i), pts(j), pts(k)};

        // get subset from starters
        std::vector<vec> subset;
        for (const int s_i: starters[idx]) {
            subset.push_back( surface[s_i] );
        }

        const vec cp_p = LSP(p, subset);
        band.push_back(idx);

        if (cp_p.has_nan()) {
            cp_pts.push_back(p);
            std::cout << "encountered nan on " << std::endl;
            std::cout << p  << std::endl;
            // currently nan is encountered only for p = (1, 0, 0)
        } else {
            cp_pts.push_back( cp_p );
        }

    }
    std::cout << "band has size: " << band.size() << std::endl;

    // construct interpolation and laplacian matrices
    std::vector<double> pts_v (pts.size() );
    for (int i = 0; i < pts.size(); i++) {
        pts_v[i] = pts(i);
    }



    const sp_mat laplacian = laplacian_3d(band, N, n, dx);
    const sp_mat E = interp_matrix_3d(pts_v, pts_v, pts_v, cp_pts, band);
}
