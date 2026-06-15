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
    
    // phi_T stores \sum A^{m-i-1} B L(i)^T at each step
    Eigen::MatrixXd phi_T = Eigen::MatrixXd::Zero(nx, N);
    Eigen::MatrixXd A_pow = Eigen::MatrixXd::Identity(nx, nx);
    
    // Q_matrix weights the state cost per prediction step.
    // For IMP mode (nx==6): output-tracking cost y = Caug*x, so Q_matrix = Q * Caug^T * Caug.
    // This naturally penalizes the harmonic displacement states (indices 2,4) via the
    // observation matrix Caug = [1,0,1,0,1,0], enabling pre-emptive IMP cancellation.
    Eigen::MatrixXd Q_matrix = Eigen::MatrixXd::Zero(nx, nx);
    if (nx == 3) {
        Q_matrix(2, 2) = Q_;
    } else if (nx == 7) {
        Q_matrix(6, 6) = Q_; // Weight the 7th state (tracking error)
    } else if (nx == 6) {
        // Output-tracking: penalize y = Caug * x
        const Eigen::RowVectorXd& Cd = plant.getCaug();
        Q_matrix = Q_ * Cd.transpose() * Cd;
    } else {
        Q_matrix(1, 1) = Q_ * 0.1; // Motor velocity damping
        Q_matrix(3, 3) = Q_ * 10.0; // Axis velocity damping
        Q_matrix(4, 4) = Q_;       // Axis position tracking
    }

    // q_target: direction in state space toward the reference output.
    // For output-tracking (nx==6): q_target = Q * Caug^T.
    // For position-tracking (others): q_target points at last state.
    Eigen::VectorXd q_target = Eigen::VectorXd::Zero(nx);
    if (nx == 6) {
        q_target = Q_ * plant.getCaug().transpose();
    } else {
        q_target(nx-1) = Q_;
    }

    Eigen::VectorXd current_L = network.getInitialState();

    for (int m = 1; m <= Np_; ++m) {
        // Recursive prediction update (Matrix-Vector form)
        phi_T = A_aug * phi_T + B_aug * current_L.transpose();

        A_pow = A_pow * A_aug;

        Omega += phi_T.transpose() * Q_matrix * phi_T + current_L * R_ * current_L.transpose();
        Psi_x += phi_T.transpose() * Q_matrix * A_pow;
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
    double delta_u_dist = -(disturbance_estimate - last_dist_);
    last_dist_ = disturbance_estimate;

    double delta_u = delta_u_mpc + delta_u_dist;

    // 3. Backlash Rejection (Punch-through)
    // Only fire on a SUSTAINED direction reversal -- the cumulative motor
    // command has been one-sided for at least kPunchSustained frames before
    // the reversal. Under sinusoidal PE the controller's current_u crosses
    // zero every half-period; firing on every crossing adds +/-backlash
    // bias each cycle and disrupts smooth tracking. Real backlash gap
    // traversal happens at step recoveries and on slow PE -- both produce
    // long stretches of one-sided current_u between reversals.
    if (backlash_ > 0.0) {
        // 30 frames at typical dt=2s = 60 seconds of sustained one-sided motor
        // command before the next reversal is credited as a real backlash-gap
        // traversal. Lower thresholds fire on every PE half-period and add
        // +/-backlash bias each cycle, disrupting smooth tracking (verified
        // experimentally on H13: threshold=10 -> 1.282", threshold=disabled
        // -> 1.185"). 30 is high enough to stay dormant under typical PE
        // (T = 30-480 s, half-period in frames = T / (2*dt) << 30 for dt >= 2)
        // and fires on step recoveries / direction changes triggered by real
        // disturbances.
        constexpr int kPunchSustained = 30;
        const double next_u = current_u_ + delta_u;
        const double new_sign = (next_u >  1e-6) ?  1.0
                              : (next_u < -1e-6) ? -1.0 : 0.0;
        const bool reversed = (last_current_u_sign_ != 0.0)
                            && (new_sign != 0.0)
                            && (new_sign != last_current_u_sign_);

        if (reversed && frames_since_reversal_ >= kPunchSustained) {
            delta_u += new_sign * backlash_;
        }

        if (reversed) frames_since_reversal_ = 0;
        else if (new_sign != 0.0) ++frames_since_reversal_;

        if (new_sign != 0.0) last_current_u_sign_ = new_sign;
    }

    // 4. Per-step rate limit.
    // The MPC gains assume the commanded delta_u takes effect immediately on
    // the measured output (rigid pure-integrator plant). When the actual mount
    // has unmodeled dynamics (compliance, stiction, large transport delay) the
    // controller can issue catastrophically large corrections trying to fix
    // what it sees. Capping the per-step demand prevents single-step blowups
    // while still allowing cumulative drift correction over many frames.
    // 5" per frame is well above typical guide pulses (sub-arcsec) but bounds
    // the worst-case excursion under plant-model mismatch.
    const double MAX_DELTA_U = 5.0;
    if (delta_u >  MAX_DELTA_U) delta_u =  MAX_DELTA_U;
    if (delta_u < -MAX_DELTA_U) delta_u = -MAX_DELTA_U;

    current_u_ += delta_u;
    last_delta_u_ = delta_u;

    return current_u_;
}

Eigen::VectorXd MPCSolver::predictNextState(const Eigen::VectorXd& x_aug, double delta_u, const TelescopePlant& plant) {
    // Prediction: x(k+1) = A_aug * x(k) + B_aug * delta_u
    // This is the "Digital Twin" prediction for the ModelMonitor
    return plant.getAaug() * x_aug + plant.getBaug() * delta_u;
}
