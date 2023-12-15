#pragma once

#include "FieldIntegration.h"

namespace SurfaceFields {
    // Our method that estimates s using eigenvalue problem
    class SpectralLocalIntegration : public LocalFieldIntegration
    {
    public:
        SpectralLocalIntegration(double sSmoothnessReg) : sreg_(sSmoothnessReg) {}

        void locallyIntegrateOneComponent(const Surface& surf, const Eigen::MatrixXd& v, Eigen::VectorXd& s);

    private:
        double sreg_;
    };
}