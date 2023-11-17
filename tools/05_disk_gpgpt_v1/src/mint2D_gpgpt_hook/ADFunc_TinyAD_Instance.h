#ifndef ADFUNC_TINYAD_INSTANCE_H
#define ADFUNC_TINYAD_INSTANCE_H

// #include "polyscope/polyscope.h"

// #include <TinyAD/Utils/LinearSolver.hh>

#include "ADFuncRunner.h"


#include <Eigen/Core>



/// This class wraps your tiny ad func in a templated way to allow for taking newton steps in sequence
template<typename N>
class ADFunc_TinyAD_Instance : public ADFuncRunner { 
            // ScalarFunction _f; public: 
            // FuncEvaluator(f) : _f(f) {} 
            // In practice: use f to evaluate the function, gradient, and hessian 

    ~ADFunc_TinyAD_Instance(){};
            
    double eval_func_local(const Eigen::VectorXd &x){
        _fun_val = _func.eval(x);
    }; 
    
    void eval_func_with_derivatives(const Eigen::VectorXd &x){
        // std::tuple<double, Eigen::VectorX<double>, Eigen::SparseMatrix<double>> = 
        auto [f, g, H_proj] = func.eval_with_derivatives(x);
        _fun_val = f;
        _grad = g; 
        _hess = H_proj;
    }; 
    void eval_func_and_proj_hess_to_psd_local(const Eigen::VectorXd &x){
        auto [f_h, g_h, H_proj_h] = func.eval_with_hessian_proj(x);
        _fun_val = f_h;
        _grad = g_h; 
        _hess = H_proj_h;

        // d = TinyAD::newton_direction(g, H_proj, solver);
        // dec = TinyAD::newton_decrement(d, g);
        
    }; 


    decltype(TinyAD::scalar_function<N>(TinyAD::range(1))) _func;

            // config flags 



}; 
            




#endif

