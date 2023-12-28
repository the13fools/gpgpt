// MyConfig.h

#ifndef MY_CONFIG_H
#define MY_CONFIG_H

#include <string>

class MyConfig {
    public:
        MyConfig();
        double w_bound;
        double w_smooth;
        double w_smooth_vector; // unused
        double w_curl;
        double w_curl_L2;
        double w_curl_L4;
        double w_s_perp; // unused 

        double w_attenuate;

        double convergence_eps;
        double identity_weight;
        double prev_energy;
        bool useProjHessian;
};

#endif // MY_CONFIG_H
