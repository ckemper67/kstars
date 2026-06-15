#pragma once
#include <Eigen/Dense>

class LaguerreNetwork {
public:
    // N is the length of the Laguerre network (prediction horizon)
    // a is the scaling factor [0, 1)
    LaguerreNetwork(int N, double a);

    // Get the initial L(0) state vector
    Eigen::VectorXd getInitialState() const;

    // Generate the next state L(k+1) = A_l * L(k)
    Eigen::VectorXd step(const Eigen::VectorXd& current_L) const;

    // Generates the L matrix where each column is L(k) for k=0...steps-1
    Eigen::MatrixXd generateTrajectory(int steps) const;

    int getN() const { return N_; }
    double getA() const { return a_; }

private:
    int N_;
    double a_;
    double beta_;

    Eigen::MatrixXd Al_;
    Eigen::VectorXd L0_;

    void initializeNetwork();
};
