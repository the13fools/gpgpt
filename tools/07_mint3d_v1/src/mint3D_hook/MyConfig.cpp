// // MyConfig.cpp

#include "MyConfig.h"
#include <iostream>

// Define default values for MyConfig members
MyConfig::MyConfig()
    : w_bound(1.e6),
      w_smooth(0.1),
      w_smooth_vector(1),
      w_curl(1e5),
      w_s_perp(0.),
      w_attenuate(0.1),
      convergence_eps(1e-12),
      identity_weight(0.1),
      prev_energy(0.0),
      useProjHessian(true) {
}

// // Implement any other MyConfig-related functionality here if needed
