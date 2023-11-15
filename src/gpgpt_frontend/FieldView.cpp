#include "FieldView.h"

// Utility function to convert Field_View enums to strings
std::string fieldViewToString(Field_View view) {
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
std::string fieldViewToFileStub(Field_View view) {
    switch (view) {
        case vec_norms: return "norms_vec";
        case delta_norms: return "norms_delta";
        case vec_dirch: return "dirch_vec";
        case moment_dirch: return "dirch_moment";
        case primal_curl_residual: return "curl_vec";
        case sym_curl_residual: return "curl_mom";
        case gui_free: return "Free";
        default: return "Unknown";
    }
}