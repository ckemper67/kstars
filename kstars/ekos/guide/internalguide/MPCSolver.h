#pragma once
#include <Eigen/Dense>
#include "TelescopePlant.h"
#include "LaguerreNetwork.h"

class MPCSolver {
public:
    // Np is the numerical prediction horizon over which we sum the cost
    MPCSolver(const TelescopePlant& plant, const LaguerreNetwork& network, double Q, double R, int Np = 100);

    void updateWeights(double Q, double R, const TelescopePlant& plant, const LaguerreNetwork& network);
    
    // Rebuilds the algebraic matrices
    void rebuildMatrices(const TelescopePlant& plant, const LaguerreNetwork& network);

    // x_aug is [ \Delta x, theta_axis ]^T (size 3 or 5)
    // setpoint is the desired theta tracking angle
    // disturbance_estimate is the estimated torque in Nm (optional)
    double computeDeltaU(const Eigen::VectorXd& x_aug, double setpoint, double disturbance_estimate = 0.0);

    // Real-time Model Prediction for the Monitor
    Eigen::VectorXd predictNextState(const Eigen::VectorXd& x_aug, double delta_u, const TelescopePlant& plant);

    double getCurrentU() const { return current_u_; }
    void setCurrentU(double u) { current_u_ = u; }
    const Eigen::RowVectorXd& getKx() const { return Kx_; }
    double getKr() const { return Kr_; }

private:
    double Q_; // Tracking error penalty
    double R_; // Control variation penalty
    int Np_;   // Simulation horizon for cost summation

    // Precomputed optimal feedback gains (Size 3 or 5)
    Eigen::RowVectorXd Kx_;
    double Kr_;

    double current_u_ = 0.0;
    double last_delta_u_ = 0.0;
    double backlash_ = 0.0;
    double last_dist_ = 0.0;
};
