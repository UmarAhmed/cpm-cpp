#ifndef PT_CLOUD_H
#define PT_CLOUD_H


#include <armadillo>
#include <vector>


/*
 * Implementation of projection algorithm
*/

std::vector<std::vector<int>> warm_start_3d(const arma::vec pts, const std::vector<arma::vec> surface, const double threshold);


arma::vec LSP(arma::vec p, std::vector<arma::vec> pts, const double epsilon = 0.0001, const int max_steps = 5);


#endif
