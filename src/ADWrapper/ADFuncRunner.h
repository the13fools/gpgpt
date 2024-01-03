#ifndef ADFUNCRUNNER_H
#define ADFUNCRUNNER_H

// #include "polyscope/polyscope.h"

#include <TinyAD/Utils/LinearSolver.hh>


#include <Eigen/Core>

class ADFuncRunner 
{ 
    public: 
    virtual double eval_func_local(const Eigen::VectorXd &x) = 0; 
    
    virtual void eval_func_with_derivatives(const Eigen::VectorXd &x) = 0; 
    virtual void eval_func_and_proj_hess_to_psd_local(const Eigen::VectorXd &x) = 0; 

    virtual int get_num_vars() = 0;


    ADFuncRunner() = default;
    virtual ~ADFuncRunner() {};

    // virtual Eigen::VectorXd eval_grad_local(const Eigen::VectorXd &x) = 0; 
    // virtual Eigen::MatrixXd eval_hess_local(const Eigen::VectorXd &x) = 0; 

    double eval_func_at(const Eigen::VectorXd &x) { return eval_func_local(x); } 
 

// can copy tinyad and return a tuple here for convenience.
    void eval_func_with_derivatives_at(const Eigen::VectorXd &x) { eval_func_with_derivatives(x);  } 
    void eval_func_and_proj_hess_to_psd_at(const Eigen::VectorXd &x) { eval_func_and_proj_hess_to_psd_local(x);  } 
    void reset_diag_hessian_reg() { identity_weight = 1e-8; }


    Eigen::VectorXd get_current_x() { return _cur_x; } 
    void set_current_x(const Eigen::VectorXd &x) { _cur_x = x; }
    double get_fval_at_x() { return _fun_val; } 
    Eigen::VectorXd get_grad_at_x() { return _grad; } 
    Eigen::SparseMatrix<double> get_hessian_at_x() { return _hess; } 


    Eigen::VectorXd take_newton_step(const Eigen::VectorXd &x);
    void reset_params();


//// Mostly shouldn't touch this but if you want to speed up convergence can do your own regualization updates
    void set_identity_weight(double set) { identity_weight = set; }
    double set_identity_weight() { return identity_weight; }

    // void clear_state()
    // {
    //     _cur_x = _cur_x * 0;
    //     _fun_val = -1;
    //     _grad.clear();
    //     _hess.clear();
    // }


    public: 
        bool useProjHessian = true;
        double prev_energy = -100000.;
        double identity_weight = 1e-8;  // This term controls how much identity we add to the hessian to make it psd.  

        // current state  
        Eigen::VectorXd _cur_x;
        Eigen::VectorXd _newton_dir;
        double _dec;
        double _max_gradient_norm; 
        double _prev_step_progress;
        double _line_search_step_size;
        double _prev_step_time;
        double _solve_residual; 
        double _rhs_norm; 
        


    protected: 


        // runner state 
        // bool x_curr_is_new = false;
        TinyAD::LinearSolver<double> solver; // make this changable 
        // Eigen::SimplicialLDLT<Eigen::SparseMatrix<PassiveT>>>
        // #include <Eigen/SparseCholesky>
        // Eigen::CholmodSupernodalLLT< Eigen::SparseMatrix<PassiveT>> >


        // cached quantities 

        double _fun_val; 
        Eigen::VectorXd _grad;
        Eigen::SparseMatrix<double> _hess;
        Eigen::SparseMatrix<double> _hess_proj;
}; 


            

///// TODO:  Add enzyme connector here  
// template <size_t N> class ADFunc_enzyme_Instance : public ADFuncRunner { 
//             // ScalarFunction _f; public: 
//             // FuncEvaluator(f) : _f(f) {} 
//             // In practice: use f to evaluate the function, gradient, and hessian 
            
//     double eval_func_local(const VectorXd &x) override { return x.squaredNorm(); } 
//     VectorXd eval_grad_local(const VectorXd &x) override { return Matrix::Ones(); } 
//     MatrixXd eval_hess_local(const VectorXd &x) override { return Matrix::Identity(x.size(), x.size()); } 


//     decltype(TinyAD::scalar_function<N>(TinyAD::range(1))) _func;


// }; 




#endif




    // double eval_func(const Eigen::VectorXd &x) { return eval_func_local(x); } 
    // Eigen::VectorXd eval_grad(const Eigen::VectorXd &x) { return eval_grad_local(x); } 
    // Eigen::MatrixXd eval_hess(const Eigen::VectorXd &x) { return eval_hess_local(x); } 
        
    // double eval_func(const Eigen::VectorXd &x) { return eval_func_local(x); } 
    // Eigen::VectorXd eval_grad(const Eigen::VectorXd &x) { return eval_grad_local(x); } 
    // Eigen::MatrixXd eval_hess(const Eigen::VectorXd &x) { return eval_hess_local(x); } 
