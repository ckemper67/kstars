#include "LaguerreNetwork.h"
#include <cmath>

LaguerreNetwork::LaguerreNetwork(int N, double a) : N_(N), a_(a) {
    beta_ = std::sqrt(1.0 - a_*a_);
    initializeNetwork();
}

void LaguerreNetwork::initializeNetwork() {
    Al_ = Eigen::MatrixXd::Zero(N_, N_);
    L0_ = Eigen::VectorXd::Zero(N_);
    
    // Initial state L(0)
    double term = beta_;
    for (int i = 0; i < N_; ++i) {
        L0_(i) = term;
        term *= -a_;
    }

    // Build the dynamic matrix A_l
    for (int i = 0; i < N_; ++i) {
        for (int j = 0; j <= i; ++j) {
            if (i == j) {
                Al_(i, j) = a_;
            } else {
                Al_(i, j) = beta_ * beta_ * std::pow(-a_, i - j - 1);
            }
        }
    }
}

Eigen::VectorXd LaguerreNetwork::getInitialState() const {
    return L0_;
}

Eigen::VectorXd LaguerreNetwork::step(const Eigen::VectorXd& current_L) const {
    return Al_ * current_L;
}

Eigen::MatrixXd LaguerreNetwork::generateTrajectory(int steps) const {
    Eigen::MatrixXd trajectory(N_, steps);
    Eigen::VectorXd current_L = L0_;
    for (int i = 0; i < steps; ++i) {
        trajectory.col(i) = current_L;
        current_L = Al_ * current_L;
    }
    return trajectory;
}
