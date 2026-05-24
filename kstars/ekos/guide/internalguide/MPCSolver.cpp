#include "MPCSolver.h"
#include "svdsolve.c"
#include <iostream>
#include <vector>
#include <iomanip>

MPCSolver::MPCSolver(const TelescopePlant& plant, const LaguerreNetwork& network, double Q, double R, int Np)
    : Q_(Q), R_(R), Np_(Np) {
    rebuildMatrices(plant, network);
}

void MPCSolver::updateWeights(double Q, double R, const TelescopePlant& plant, const LaguerreNetwork& network) {
    Q_ = Q;
    R_ = R;
    rebuildMatrices(plant, network);
}

void MPCSolver::rebuildMatrices(const TelescopePlant& plant, const LaguerreNetwork& network) {
    int N = network.getN();
    int nx = plant.getOrder();
    
    backlash_ = plant.getBacklash();

    Eigen::MatrixXd Omega = Eigen::MatrixXd::Zero(N, N);
    Eigen::VectorXd Psi_r = Eigen::VectorXd::Zero(N);
    Eigen::MatrixXd Psi_x = Eigen::MatrixXd::Zero(N, nx);
    
    // Get adaptive matrices from the universal plant
    const Eigen::MatrixXd& A_aug = plant.getAaug();
    const Eigen::VectorXd& B_aug = plant.getBaug();
    const Eigen::RowVectorXd& C_aug = plant.getCaug();
    
    // phi_T stores \sum A^{m-i-1} B L(i)^T at each step
    Eigen::MatrixXd phi_T = Eigen::MatrixXd::Zero(nx, N);
    Eigen::MatrixXd A_pow = Eigen::MatrixXd::Identity(nx, nx);
    
    // Define a state weighting matrix Q_matrix
    // For Rigid N=3: [0, 0, Q] - only position
    // For Flexible N=5: [0, 0.1, 0, 0.5, Q] - dampen velocities
    Eigen::MatrixXd Q_matrix = Eigen::MatrixXd::Zero(nx, nx);
    if (nx == 3) {
        Q_matrix(2, 2) = Q_;
    } else {
        Q_matrix(1, 1) = Q_ * 0.1; // Motor velocity damping
        Q_matrix(3, 3) = Q_ * 10.0; // Axis velocity damping (Aggressive vibration suppression)
        Q_matrix(4, 4) = Q_;       // Axis position tracking
    }
    
    Eigen::VectorXd current_L = network.getInitialState();
    
    for (int m = 1; m <= Np_; ++m) {
        // Recursive prediction update (Matrix-Vector form)
        phi_T = A_aug * phi_T + B_aug * current_L.transpose();
        
        A_pow = A_pow * A_aug;
        
        // Cost = Predict(x)^T * Q * Predict(x)
        // Predict(x) = A_pow * x + phi_T * coefficients
        // Contribution to Omega (Hessian): phi_T^T * Q * phi_T
        Omega += phi_T.transpose() * Q_matrix * phi_T + current_L * R_ * current_L.transpose();
        
        // Contribution to Psi_x (State feedback): phi_T^T * Q * A_pow
        Psi_x += phi_T.transpose() * Q_matrix * A_pow;

        // Contribution to Psi_r (Reference tracking): phi_T^T * Q_vec * target
        // Since setpoint only applies to position (last state)
        Eigen::VectorXd q_target = Eigen::VectorXd::Zero(nx);
        q_target(nx-1) = Q_;
        Psi_r += phi_T.transpose() * q_target;

        current_L = network.step(current_L);
    }
    
    // Solve for the optimal gain matrix K = inv(Omega) * Psi
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Omega_rm = Omega;
    std::vector<double> w(N);
    std::vector<double> v(N * N);
    
    // SVD decomposition (Embedded-safe, zero allocation)
    svd(N, N, Omega_rm.data(), w.data(), v.data());
    
    // Solve for reference gain Kr
    Eigen::VectorXd K_r_vec(N);
    svd_backsub(N, N, Omega_rm.data(), w.data(), v.data(), Psi_r.data(), 1e-12, K_r_vec.data());
    
    // Solve for state gains Kx (column-by-column)
    Eigen::MatrixXd K_xT(N, nx);
    for (int i = 0; i < nx; ++i) {
        Eigen::VectorXd b_col = Psi_x.col(i);
        Eigen::VectorXd x_col(N);
        svd_backsub(N, N, Omega_rm.data(), w.data(), v.data(), b_col.data(), 1e-12, x_col.data());
        K_xT.col(i) = x_col;
    }
    
    // Collapse Laguerre coefficients into the hot-loop Gains
    Eigen::VectorXd L0 = network.getInitialState();
    Kr_ = L0.dot(K_r_vec);
    Kx_ = L0.transpose() * K_xT;
}

double MPCSolver::computeDeltaU(const Eigen::VectorXd& x_aug, double setpoint, double disturbance_estimate) {
    // 1. Calculate dynamic MPC correction
    double delta_u_mpc = Kr_ * setpoint - Kx_.dot(x_aug);
    
    // 2. Disturbance Feedforward (DOB)
    static double last_dist = 0.0;
    double delta_u_dist = -(disturbance_estimate - last_dist);
    last_dist = disturbance_estimate;

    double delta_u = delta_u_mpc + delta_u_dist;

    // 3. Backlash Rejection (Punch-through)
    // If the demand reverses the sign of the current effort, we must cross the deadzone
    if (backlash_ > 0.0 && std::abs(current_u_) > 1e-6) {
        bool dir_changed = (current_u_ > 0 && (current_u_ + delta_u) < 0) ||
                           (current_u_ < 0 && (current_u_ + delta_u) > 0);
        if (dir_changed) {
            // Add a one-time bias to hop over the gear play
            delta_u += (delta_u > 0 ? 1.0 : -1.0) * backlash_;
        }
    }

    current_u_ += delta_u;
    last_delta_u_ = delta_u;
    
    return current_u_;
}

Eigen::VectorXd MPCSolver::predictNextState(const Eigen::VectorXd& x_aug, double delta_u, const TelescopePlant& plant) {
    // Prediction: x(k+1) = A_aug * x(k) + B_aug * delta_u
    // This is the "Digital Twin" prediction for the ModelMonitor
    return plant.getAaug() * x_aug + plant.getBaug() * delta_u;
}
