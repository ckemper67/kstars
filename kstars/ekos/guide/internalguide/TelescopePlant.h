#pragma once
#include <Eigen/Dense>

class TelescopePlant {
public:
    TelescopePlant(double J_motor, double Bf, double Kt, double dt);

    void updateParameters(double J_motor, double Bf, double Kt, double dt);
    
    // Universal Mode setters
    void setAxisInertia(double J_axis) { J_axis_ = J_axis; discretize(); }
    void setSpringStiffness(double Ks) { Ks_ = Ks; discretize(); }
    void setSpringDamping(double Bs) { Bs_ = Bs; discretize(); }
    void setBacklash(double delta) { backlash_ = delta; }

    // Adaptive matrix access (size depends on getOrder())
    const Eigen::MatrixXd& getAaug() const { return Aaug_; }
    const Eigen::VectorXd& getBaug() const { return Baug_; }
    const Eigen::RowVectorXd& getCaug() const { return Caug_; }
    
    // Raw Motor matrices (2D)
    const Eigen::Matrix2d& getAd() const { return Ad_; }
    const Eigen::Vector2d& getBd() const { return Bd_; }

    int getOrder() const { return system_order_; }
    double getDt() const { return dt_; }
    double getBacklash() const { return backlash_; }
    
    void setCoulombFriction(double Tc) { Tc_ = Tc; }
    double getCoulombFriction() const { return Tc_; }

private:
    double J_motor_;
    double J_axis_ = 0.0;
    double Bf_;
    double Kt_;
    double dt_;

    double Ks_ = 1e9; // Default to rigid (High stiffness)
    double Bs_ = 0.1;
    double Tc_ = 0.0;
    double backlash_ = 0.0;

    int system_order_ = 3;

    // Base motor matrices
    Eigen::Matrix2d Ad_;
    Eigen::Vector2d Bd_;

    // Augmented MPC matrices (Adaptive size)
    Eigen::MatrixXd Aaug_;
    Eigen::VectorXd Baug_;
    Eigen::RowVectorXd Caug_;

    void discretize();
};
