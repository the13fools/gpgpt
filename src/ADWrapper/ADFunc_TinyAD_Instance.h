#ifndef ADFUNC_TINYAD_INSTANCE_H
#define ADFUNC_TINYAD_INSTANCE_H

// #include "polyscope/polyscope.h"

// #include <TinyAD/Utils/LinearSolver.hh>

#include "ADWrapper/ADFuncRunner.h"
#include <Eigen/Core>

#include <TinyAD/ScalarFunction.hh>
// #include <TinyAD/Utils/LinearSolver.hh>


/// This class wraps your tiny ad func in a templated way to allow for taking newton steps in sequence
template<int N>
class ADFunc_TinyAD_Instance : public ADFuncRunner { 
            // ScalarFunction _f; public: 
            // FuncEvaluator(f) : _f(f) {} 
            // In practice: use f to evaluate the function, gradient, and hessian 
    public: 

    ADFunc_TinyAD_Instance(){
        _cur_x = Eigen::VectorXd::Random(N) * 1e-5;
    };
    ~ADFunc_TinyAD_Instance(){};

    // ADFunc_TinyAD_Instance& operator=(const ADFunc_TinyAD_Instance& other)
    // {
    //     if (this != &other) 
    //     {
    //         this->set_tinyad_objective_func(other._func);
    //         // Copy each member from 'other' to 'this'
    //         // Be careful with pointers! You might need to perform a deep copy
    //     }
    //     return *this;
    // }
                
    double eval_func_local(const Eigen::VectorXd &x){
        _fun_val = _func->eval(x);
        return _fun_val;
    }; 
    
    void eval_func_with_derivatives(const Eigen::VectorXd &x){
        // std::tuple<double, Eigen::VectorX<double>, Eigen::SparseMatrix<double>> = 
        auto [f, g, H_proj] = _func->eval_with_derivatives(x);
        _fun_val = f;
        _grad = g; 
        _hess = H_proj;

        // _max_gradient_norm = g.cwiseAbs().maxCoeff(); 
    }; 
    void eval_func_and_proj_hess_to_psd_local(const Eigen::VectorXd &x){
        auto [f_h, g_h, H_proj_h] = _func->eval_with_hessian_proj(x);
        _fun_val = f_h;
        _grad = g_h; 
        _hess = H_proj_h;

        // _max_gradient_norm = g.cwiseAbs().maxCoeff(); 

        // d = TinyAD::newton_direction(g, H_proj, solver);
        // dec = TinyAD::newton_decrement(d, g);
        
    }; 

    void set_tinyad_objective_func(decltype(TinyAD::scalar_function<N>(TinyAD::range(1)))* func)
    {
        _func = func;
        _cur_x = Eigen::VectorXd::Zero(_func->n_vars);
    }

    int get_num_vars()
    {
        // std::cout << std::endl << "get_num_vars() called" << std::endl;
        // std::cout << "_func->n_vars = " << _func->n_vars << std::endl;
        return _func->n_vars;
    }


    decltype(TinyAD::scalar_function<N>(TinyAD::range(1)))* _func;

            // config flags 



}; 
            




#endif



        // ADFunc_TinyAD_Instance(decltype(TinyAD::scalar_function<N>(TinyAD::range(1))) func) : _func(func) {
            // std::cout << "ADFunc_TinyAD_Instance constructor called" << std::endl;
            // std::cout << "func = " << func << std::endl;
            // std::cout << "_func = " << _func << std::endl;
            // std::cout << "_func.eval(Eigen::VectorXd::Zero(N)) = " << _func.eval(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_derivatives(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_derivatives(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_hessian_proj(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_hessian_proj(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_hessian(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_hessian(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_hessian(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_hessian(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_hessian(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_hessian(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_hessian(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_hessian(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_hessian(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_hessian(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_hessian(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_hessian(Eigen::VectorXd::Zero(N)) << std::endl;
            // std::cout << "_func.eval_with_hessian(Eigen::VectorXd::Zero(N)) = " << _func.eval_with_hessian(Eigen::VectorXd::Zero(N
        // };