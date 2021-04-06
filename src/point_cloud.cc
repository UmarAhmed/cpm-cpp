#include <iostream>
#include <armadillo>
#include <vector>
#include <algorithm>


using namespace std;
using namespace arma;

/*
 * vec = [x0, ..., xn] of the cartesian product grid
*/
std::vector<std::vector<int>> warm_start_3d(const vec pts, const std::vector<vec> surface, const double threshold) {
    const int n = pts.size();
    std::vector<std::vector<int>> grid_dict(n * n * n);
    const double dx = pts(1) - pts(0);
    const int r = ceil(threshold / dx);

    for (int i = 0; i < surface.size(); i++) {
        const vec s = surface[i];

        const int grid_x = round( (s(0) - pts(0)) / dx);
        const int grid_y = round( (s(1) - pts(0)) / dx);
        const int grid_z = round( (s(2) - pts(0)) / dx);

        const int leftx = std::max(0, grid_x - r);
        const int rightx = std::min(n - 1, grid_x + r);


        const int lefty = std::max(0, grid_y - r);
        const int righty = std::min(n - 1, grid_y + r);

        const int leftz = std::max(0, grid_z - r);
        const int rightz = std::min(n - 1, grid_z + r);

        for (int row = leftx; row < rightx + 1; row++) {
            for (int col = lefty; col < righty + 1; col++) {
                for (int za = leftz; za < rightz + 1; za++) {
                    const vec q {pts(row), pts(col), pts(za)};
                    if (norm(q - s) <= threshold) {
                        const int idx = row + n * (col + n * za);
                        grid_dict[idx].emplace_back(i);
                    }
                    // only add points in band
                }
            }
        }
    }
    return grid_dict;
}

/*
 * Implementation of Automatic-Least Square Projection
*/
vec LSP(vec p, std::vector<vec> pts, const double epsilon = 0.0001, const int max_steps = 5) {
    int k = 0;
    double t = 0;
    vec n;
    
    while (k < max_steps) {
        // Compute weights
        vec weights = vec(pts.size(), fill::zeros);
        for (int i = 0; i < pts.size(); i++) {
            weights(i) = 1 / (1 + pow(norm(p - pts[i]), 4) );
        }

        // Compute projection direction
        int dim = p.size();
        double c_0 = sum(weights);
        vec c = vec(dim, fill::zeros);
        for (int d = 0; d < dim; d++) {
            for (int i = 0; i < pts.size(); i++) {
                c(d) += weights[i] * pts[i](d);
            }
        }
        vec m = (c / c_0) - p;
        n = m / norm(m);

        // Compute projection through Directed Projection
        double tnew = dot(c, n) / c_0 - dot(p, n);

        // Check and update t
        if ( abs(tnew - t) < epsilon) {
            break;
        }
        t = tnew;


        // Update set of points by looking at weights
        double w_max = weights.max();
        double w_mean = c_0 / weights.size();
        double w_lim = w_mean;
        if (k < 11) {
            w_lim += (w_max - w_mean) / (12 - k);
        } else {
            w_lim += (w_max - w_mean);
        }

        vector<vec> newpts;
        for (int i = 0; i < pts.size(); i++) {
            if (weights(i) >= w_lim) {
                newpts.emplace_back(std::move(pts[i]));
            }
        }

        if (newpts.size() == 0) {
            break;
        }

        pts = newpts;
        k++;
    }
    return p + t * n;
}

/*
int main() {

    vector<vec> pts;
    double t = 0;
    while (t <= 6.28) {
        vec pt = {sin(t), cos(t)};
        pts.push_back(pt);
        t += 0.000628;
    }


    wall_clock timer;
    timer.tic();

    auto N = 10000;
    for (int i = 0; i < N; i++) {
        vec p(2, fill::randu);
        auto result = LSP(p, pts);
    }

    double n = timer.toc();
    cout << "Projecting " << N << " points onto a set of " << pts.size() << " took " << n << endl; 
}
*/
