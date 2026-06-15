/*
    SPDX-FileCopyrightText: 2026 Christian Kemper <ckemper@gmail.com>

    SPDX-License-Identifier: GPL-2.0-or-later
*/

#include "fftperiodestimator.h"
#include "ekos/guide/internalguide/MPI_IS_gaussian_process/src/math_tools.h"
#include <cmath>
#include <algorithm>
#include <numeric>

FFTPeriodEstimator::FFTPeriodEstimator(int maxBufferSize, double minPeriod, double maxPeriod)
    : maxBufferSize_(maxBufferSize), minPeriod_(minPeriod), maxPeriod_(maxPeriod) {
}

void FFTPeriodEstimator::addDataPoint(double timestamp, double openLoopError) {
    timestamps_.push_back(timestamp);
    errors_.push_back(openLoopError);

    if (static_cast<int>(timestamps_.size()) > maxBufferSize_) {
        timestamps_.erase(timestamps_.begin());
        errors_.erase(errors_.begin());
    }
}

void FFTPeriodEstimator::reset() {
    timestamps_.clear();
    errors_.clear();
}

bool FFTPeriodEstimator::estimatePeriods(double currentPeriod, double& learnedT1, double& learnedT2) {
    int n = static_cast<int>(timestamps_.size());
    if (n < 40) {
        return false;
    }

    double totalTime = timestamps_.back() - timestamps_.front();
    if (totalTime < 2.0 * currentPeriod) {
        return false;
    }

    Eigen::VectorXd time(n);
    Eigen::VectorXd data(n);
    for (int i = 0; i < n; ++i) {
        time(i) = timestamps_[i];
        data(i) = errors_[i];
    }

    // 1. Linear detrending
    double mean_x = time.mean();
    double mean_y = data.mean();
    double num = 0.0;
    double den = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = time(i) - mean_x;
        num += dx * (data(i) - mean_y);
        den += dx * dx;
    }
    double slope = (den > 1e-9) ? (num / den) : 0.0;
    double intercept = mean_y - slope * mean_x;
    Eigen::VectorXd detrended_data = data - (slope * time + Eigen::VectorXd::Constant(n, intercept));

    // 2. Hamming window
    Eigen::VectorXd windowed_data = detrended_data.array() * math_tools::hamming_window(n).array();

    // 3. Compute spectrum using zero-padding size of 4096
    std::pair<Eigen::VectorXd, Eigen::VectorXd> result = math_tools::compute_spectrum(windowed_data, 4096);
    if (result.first.size() == 0 || result.second.size() == 0) {
        return false;
    }

    Eigen::ArrayXd amplitudes = result.first;
    Eigen::ArrayXd frequencies = result.second;

    // 4. Correct for average time step width
    double dt = totalTime / (n - 1);
    frequencies /= dt;

    // 5. Filter out invalid periods
    Eigen::ArrayXd periods = 1.0 / frequencies;
    for (int i = 0; i < amplitudes.size(); ++i) {
        double p = periods(i);
        if (p < minPeriod_ || p > maxPeriod_) {
            amplitudes(i) = 0.0;
        }
    }

    // 6. Find primary peak (fundamental period)
    Eigen::VectorXd::Index maxIndex1;
    amplitudes.maxCoeff(&maxIndex1);
    if (amplitudes(maxIndex1) < 1e-5) {
        return false;
    }

    double freq1 = performQuadraticInterpolation(amplitudes, frequencies, maxIndex1);
    if (freq1 <= 0.0) {
        return false;
    }
    learnedT1 = 1.0 / freq1;

    // 7. Find secondary peak (harmonic period)
    Eigen::ArrayXd amplitudes2 = amplitudes;
    int win = 50; // zero-out window around primary peak (prevents spectral leakage cross-talk)
    int start = std::max(0, static_cast<int>(maxIndex1) - win);
    int end = std::min(static_cast<int>(amplitudes2.size()) - 1, static_cast<int>(maxIndex1) + win);
    for (int i = start; i <= end; ++i) {
        amplitudes2(i) = 0.0;
    }

    Eigen::VectorXd::Index maxIndex2;
    amplitudes2.maxCoeff(&maxIndex2);
    if (amplitudes2(maxIndex2) > 0.15 * amplitudes(maxIndex1)) {
        double freq2 = performQuadraticInterpolation(amplitudes2, frequencies, maxIndex2);
        if (freq2 > 0.0) {
            learnedT2 = 1.0 / freq2;
        } else {
            learnedT2 = learnedT1 / 2.0;
        }
    } else {
        learnedT2 = 0.0; // No significant secondary harmonic
    }

    return true;
}

double FFTPeriodEstimator::performQuadraticInterpolation(const Eigen::ArrayXd& amplitudes, const Eigen::ArrayXd& frequencies, int maxIndex) {
    double max_frequency = frequencies(maxIndex);

    if (maxIndex < frequencies.size() - 1 && maxIndex > 0) {
        double spread = std::abs(frequencies(maxIndex - 1) - frequencies(maxIndex + 1));

        Eigen::VectorXd interp_loc(3);
        interp_loc << frequencies(maxIndex - 1), frequencies(maxIndex), frequencies(maxIndex + 1);
        interp_loc = interp_loc.array() - max_frequency; // centering for numerical stability
        interp_loc = interp_loc.array() / spread;        // normalize

        Eigen::VectorXd interp_dat(3);
        interp_dat << amplitudes(maxIndex - 1), amplitudes(maxIndex), amplitudes(maxIndex + 1);
        interp_dat = interp_dat.array() / amplitudes(maxIndex); // normalize

        if (interp_dat.maxCoeff() - interp_dat.minCoeff() >= 1e-10) {
            // Build quadratic fitting matrix
            Eigen::MatrixXd phi(3, 3);
            phi.row(0) = interp_loc.array().pow(2);
            phi.row(1) = interp_loc.array().pow(1);
            phi.row(2) = interp_loc.array().pow(0);

            // Linear regression
            Eigen::VectorXd w = (phi * phi.transpose()).ldlt().solve(phi * interp_dat);

            // Recover interpolated frequency peak
            max_frequency = max_frequency - w(1) / (2.0 * w(0)) * spread;
        }
    }

    return max_frequency;
}
