#ifndef GNGLOBALINTEGRATION_H
#define GNGLOBALINTEGRATION_H

#include "GlobalFieldIntegration.h"

namespace SurfaceFields {
class MIGlobalIntegration : public  GlobalFieldIntegration
{
public:
    MIGlobalIntegration(double anisotropy, double smoothnessReg) : aniso_(anisotropy),
        smoothreg_(smoothnessReg) {}

    virtual void globallyIntegrateOneComponent(const Surface& surf, const Eigen::MatrixXd& v, Eigen::VectorXd& theta, Eigen::VectorXd* edge_omega = nullptr) override;
    
private:
    double aniso_;
    double smoothreg_;
};
}

#endif
