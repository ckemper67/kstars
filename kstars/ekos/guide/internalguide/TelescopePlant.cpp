#include "TelescopePlant.h"
#include <cmath>
#include <Eigen/Dense>

TelescopePlant::TelescopePlant(double J_motor, double Bf, double Kt, double dt)
    : J_motor_(J_motor), Bf_(Bf), Kt_(Kt), dt_(dt) {
    discretize();
}

void TelescopePlant::updateParameters(double J_motor, double Bf, double Kt, double dt) {
    J_motor_ = J_motor;
    Bf_ = Bf;
    Kt_ = Kt;
    dt_ = dt;
    discretize();
}

void TelescopePlant::discretize() {
    Ad_ = Eigen::Matrix2d::Identity();
    Bd_ = Eigen::Vector2d::Zero();

    // Determine if we should use the Rigid (3-state) or Flexible (5-state) model
    if (Ks_ > 1e7 || J_axis_ < 1e-6) {
        // --- RIGID MODE (3-state MPC) ---
        system_order_ = 3;
        
        double p = Bf_ / J_motor_;
        if (std::abs(p) < 1e-6) {
            Ad_ << 1.0, dt_,
                   0.0, 1.0;
            Bd_ << 0.5 * Kt_ * dt_ * dt_ / J_motor_,
                   Kt_ * dt_ / J_motor_;
        } else {
            double exp_pt = std::exp(-p * dt_);
            Ad_ << 1.0, (1.0 - exp_pt) / p,
                   0.0, exp_pt;
            Bd_ << Kt_ * (dt_ / p - (1.0 - exp_pt) / (p * p)) / J_motor_,
                   Kt_ * (1.0 - exp_pt) / (p * J_motor_);
        }

        Aaug_ = Eigen::MatrixXd::Zero(3, 3);
        Aaug_.block<2,2>(0,0) = Ad_;
        Aaug_.block<1,2>(2,0) = Ad_.row(0); // theta(k+1) is the first row of Ad * x
        Aaug_(2,2) = 1.0;

        Baug_ = Eigen::VectorXd::Zero(3);
        Baug_.segment<2>(0) = Bd_;
        Baug_(2) = Bd_(0);

        Caug_ = Eigen::RowVectorXd::Zero(3);
        Caug_(2) = 1.0;
    } else {
        // --- FLEXIBLE MODE (5-state MPC) ---
        system_order_ = 5;

        // Continuous matrices for 2-mass system: [theta_m, v_m, theta_a, v_a]
        Eigen::Matrix4d Ac = Eigen::Matrix4d::Zero();
        Ac(0, 1) = 1.0;
        Ac(1, 0) = -Ks_ / J_motor_;
        Ac(1, 1) = -(Bf_ + Bs_) / J_motor_;
        Ac(1, 2) = Ks_ / J_motor_;
        Ac(1, 3) = Bs_ / J_motor_;
        Ac(2, 3) = 1.0;
        Ac(3, 0) = Ks_ / J_axis_;
        Ac(3, 1) = Bs_ / J_axis_;
        Ac(3, 2) = -Ks_ / J_axis_;
        Ac(3, 3) = -(Bs_) / J_axis_;

        Eigen::Vector4d Bc = Eigen::Vector4d::Zero();
        Bc(1) = Kt_ / J_motor_;

        // Discretize using 2nd-order Taylor expansion
        Eigen::Matrix4d I = Eigen::Matrix4d::Identity();
        Eigen::Matrix4d Ad_4 = I + Ac * dt_ + 0.5 * Ac * Ac * dt_ * dt_;
        Eigen::Vector4d Bd_4 = (I * dt_ + 0.5 * Ac * dt_ * dt_) * Bc;

        // Save motor-side 2D matrices for observer use
        Ad_ = Ad_4.block<2,2>(0,0);
        Bd_ = Bd_4.segment<2>(0);

        // Build 5D Augmented System: [del_theta_m, del_v_m, del_theta_a, del_v_a, theta_a]
        Aaug_ = Eigen::MatrixXd::Zero(5, 5);
        Aaug_.block<4,4>(0,0) = Ad_4;
        
        // Update Axis Position (theta_a is the 3rd state, index 2, in Ad_4)
        Aaug_.block<1,4>(4,0) = Ad_4.row(2);
        Aaug_(4,4) = 1.0;

        Baug_ = Eigen::VectorXd::Zero(5);
        Baug_.segment<4>(0) = Bd_4;
        Baug_(4) = Bd_4(2);

        Caug_ = Eigen::RowVectorXd::Zero(5);
        Caug_(4) = 1.0;
        /*
        std::cout << "[UniversalPlant] Flexible 5-state matrices built.\n";
        std::cout << "Ad_4 (Motor+Axis):\n" << Ad_4 << "\n";
        std::cout << "Bd_4 (Motor+Axis):\n" << Bd_4.transpose() << "\n";
        */
    }
}
