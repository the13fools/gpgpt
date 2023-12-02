#include "FieldView.h"

// Utility function to convert Field_View enums to strings
std::string Views::fieldViewToString(Field_View view) {
    switch (view) {
        case vec_norms: return "Vector Norms";
        case delta_norms: return "Delta Norms";
        case vec_dirch: return "Vec Dirichlet";
        case moment_dirch: return "Moment Dirichlet";
        case primal_curl_residual: return "Primal Curl Residual";
        case sym_curl_residual: return "Symmetric Curl Residual";
        case gui_free: return "Free";
        default: return "Unknown";
    }
}

// Utility function to convert Field_View enums to strings
std::string Views::fieldViewToFileStub(Field_View view) {
    switch (view) {
        case vec_norms: return "nv_norms_vec";
        case delta_norms: return "nd_norms_delta";
        case vec_dirch: return "dv_dirch_vec";
        case moment_dirch: return "dm_dirch_moment";
        case primal_curl_residual: return "cv_curl_vec";
        case sym_curl_residual: return "cm_curl_mom";
        case gui_free: return "free";
        default: return "Unknown";
    }
}