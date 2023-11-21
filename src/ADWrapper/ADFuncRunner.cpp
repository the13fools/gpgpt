    
    
    
    #include "ADFuncRunner.h"
    
    #include <TinyAD/Utils/LinearSolver.hh>

    #include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

// #include <fstream>
#include <sys/stat.h>
#include <iostream>
    
    
    void ADFuncRunner::reset_params()
    {
        useProjHessian = true;
        prev_energy = -100000.;
        identity_weight = 1e-8; 
    }
    
    
    Eigen::VectorXd ADFuncRunner::take_newton_step(const Eigen::VectorXd &x){
                   


            // auto [f, g, H_proj] = func.eval_with_hessian_proj(x);
            // auto [f, g, H_proj] = func.eval_with_derivatives(x);

            this->eval_func_with_derivatives(x);
            double f = this->get_fval_at_x();
            Eigen::VectorXd g = this->get_grad_at_x();
            Eigen::SparseMatrix<double> H = this->get_hessian_at_x();

            Eigen::SparseMatrix<double> H_proj;

            // std::cout<<"the number of nonzeros "<<H_proj.nonZeros() << "number of non-zeros per dof " << H_proj.nonZeros() / (6*F.rows()) << " # rows " << H_proj.rows() << " faces " << F.rows() <<std::endl;

            // std::cout<<"the number of nonzeros "<<H_proj.nonZeros()<<std::endl;
            Eigen::VectorXd d;
            double dec;
            // d = TinyAD::newton_direction(g, H_proj, solver);
             // = TinyAD::newton_decrement(d, g);

            if (prev_energy < 0)
            {
              prev_energy = f * 1.00000000001;
            }
            
            try
            {
              // First try to converge with projected hessian 
              if (useProjHessian) /// w_smooth_vector > 0 || 
              {
                this->eval_func_and_proj_hess_to_psd_local(x);
                f = this->get_fval_at_x(); // maybe don't need these two lines 
                g = this->get_grad_at_x();
                H_proj = this->get_hessian_at_x();
                d = TinyAD::newton_direction(g, H_proj, solver, 0.);
                dec = TinyAD::newton_decrement(d, g);


                // Decide when to switch to true hessian 
                if ( dec / f < 1e-3)
                {
                  useProjHessian = false;
                  std::cout << "*** switch off projected hessian to fine-tune result ***" << std::endl;
                }


              }
              else
              {
                d = TinyAD::newton_direction(g, H, solver, identity_weight);
                dec = TinyAD::newton_decrement(d, g);
                identity_weight = identity_weight / 2.;
              }
              
            }
            catch(const std::exception& e)
            {
               std::cout << "*** diagonally regularized hessian not PSD. failed to produce a descent direction.  Falling back on SVD projected hessian + increasing regularization weight. ***" << std::endl;

              this->eval_func_and_proj_hess_to_psd_local(x);
              f = this->get_fval_at_x(); // maybe don't need these two lines 
              g = this->get_grad_at_x();
              H_proj = this->get_hessian_at_x();
              d = TinyAD::newton_direction(g, H_proj, solver);
              dec = TinyAD::newton_decrement(d, g);
              if ( !useProjHessian )
                identity_weight = identity_weight * 10.;
              else 
                assert(false); // this should never happen
            }

            _cur_x = TinyAD::line_search(x, d, f, g, [&] (Eigen::VectorXd& point) -> double { return this->eval_func_at(point); }, 1., .8, 512, 1e-3);
            _newton_dir = d;
            // _dec = dec;
            _dec = f - this->eval_func_at(_cur_x);  // use this to track the true decrement 
            
            // 
            std::cout << "current decrement: " << _dec << " newton dec: " << dec << std::endl;

            return _cur_x;
    };