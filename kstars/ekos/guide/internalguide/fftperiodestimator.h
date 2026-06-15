/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#pragma once

#include <vector>
#include <Eigen/Dense>

class FFTPeriodEstimator {
public:
    FFTPeriodEstimator(int maxBufferSize = 600, double minPeriod = 10.0, double maxPeriod = 1000.0);

    // Adds a data point to the history
    void addDataPoint(double timestamp, double openLoopError);

    // Clears the history
    void reset();

    // Runs period estimation. Returns true if periods have updated/converged.
    bool estimatePeriods(double currentPeriod, double& learnedT1, double& learnedT2);

    int getNumPoints() const { return static_cast<int>(timestamps_.size()); }

private:
    int maxBufferSize_;
    double minPeriod_;
    double maxPeriod_;

    std::vector<double> timestamps_;
    std::vector<double> errors_;

    double performQuadraticInterpolation(const Eigen::ArrayXd& amplitudes, const Eigen::ArrayXd& frequencies, int maxIndex);
};
