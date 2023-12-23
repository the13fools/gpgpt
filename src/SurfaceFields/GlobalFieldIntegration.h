#ifndef FIELDINTEGRATION_H
#define FIELDINTEGRATION_H

#include "../Surface.h"
#include <Eigen/Core>

namespace SurfaceFields {
class GlobalFieldIntegration
{
public:
    GlobalFieldIntegration() {}
    virtual ~GlobalFieldIntegration() = default;

public:
    // Integrates a given vector field, assuming:
    // - the surface surf has one connected component
    // - the vector field v has no singularities
    // Result is a periodic function (values in [0, 2pi)) on the vertices of surf.
    virtual void globallyIntegrateOneComponent(const Surface& surf, const Eigen::MatrixXd& v, Eigen::VectorXd& theta, Eigen::VectorXd* edge_omega = nullptr) = 0;
};
}

#endif