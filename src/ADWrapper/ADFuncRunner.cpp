    
    
    
    #include "ADFuncRunner.h"
    
    #include <TinyAD/Utils/LinearSolver.hh>

    #include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

// #include <fstream>
#include <sys/stat.h>
#include <iostream>

#include <chrono>
    
    
    void ADFuncRunner::reset_params()
    {
        useProjHessian = true;
        prev_energy = -100000.;
        identity_weight = 1e-8; 
        identity_vanished = false;
    }
    
    
    Eigen::VectorXd ADFuncRunner::take_newton_step(const Eigen::VectorXd &x){
                   


            // auto [f, g, H_proj] = func.eval_with_hessian_proj(x);
            // auto [f, g, H_proj] = func.eval_with_derivatives(x);
            
            double cur_obj = this->eval_func_at(x);

            std::cout << std::endl << "take a newton step.  Cur Objective: " << cur_obj << std::endl;
            auto t1 = std::chrono::high_resolution_clock::now();



            if (cur_obj < 1e-15)
            {
              std::cout << "Tiny objective: " << cur_obj << "  Exiting. This shouldn't have been called in the first place." << std::endl;
              return x;
            }

            this->eval_func_with_derivatives(x);
            double f = this->get_fval_at_x();
            Eigen::VectorXd g = this->get_grad_at_x();
            _rhs_norm = g.squaredNorm();

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

            

            _cur_x = x; 

            
            if (useProjHessian) /// w_smooth_vector > 0 || 
            {
              std::cout << "Compute Hessian" << std::endl;
              this->eval_func_and_proj_hess_to_psd_local(x);
              f = this->get_fval_at_x(); // maybe don't need these two lines 
              g = this->get_grad_at_x();
              H_proj = this->get_hessian_at_x();
              std::cout << "Start Newton Solve, add timing here" << std::endl;
              d = TinyAD::newton_direction(g, H_proj, solver, 0.);
              dec = TinyAD::newton_decrement(d, g);

              _solve_residual = (H_proj * d + g).squaredNorm();
                
              


              // Decide when to switch to true hessian 
              if ( dec / f < 1e-4)
              {
                useProjHessian = false;
                std::cout << "*** switch off projected hessian to fine-tune result ***" << std::endl;
                std::cout << "*** *** *** newton decrement was: " << dec << " | fval was: " << f << " | ratio was: " << dec / f << "*** *** ***" << std::endl;
              }

              _cur_x = TinyAD::line_search(x, d, f, g, [&] (Eigen::VectorXd& point) -> double { return this->eval_func_at(point); }, 1., .8, 512, 1e-3);


            }
            else
            {



              std::cout << "current identity regularizer" << identity_weight << std::endl;
              try
              {

                d = TinyAD::newton_direction(g, H, solver, identity_weight);
                Eigen::VectorXd tmp = x + d;
                double obj_after_step = this->eval_func_at(tmp);
                dec = cur_obj - obj_after_step;
                _solve_residual = (H * d + g).squaredNorm();

                if ( dec > 0 )
                {
                  _cur_x = tmp; 
                  identity_weight = std::max(identity_weight / 2., identity_min);
                  if (identity_weight == identity_min)
                  {
                    identity_vanished = true;
                  }
                  else 
                  {
                    identity_vanished = false;
                  }
                }
                else 
                {
                  identity_weight *= 10.;
                }
                // dec = TinyAD::newton_decrement(d, g);
              }
              catch(const std::exception& e)
              {
                identity_weight *= 10.;
                // std::cerr << e.what() << '\n';
              }




            }





            


            
            /*
            try
            {
              // First try to converge with projected hessian 
              if (useProjHessian) /// w_smooth_vector > 0 || 
              {
                std::cout << "Compute Hessian" << std::endl;
                this->eval_func_and_proj_hess_to_psd_local(x);
                f = this->get_fval_at_x(); // maybe don't need these two lines 
                g = this->get_grad_at_x();
                H_proj = this->get_hessian_at_x();
                std::cout << "Start Newton Solve, add timing here" << std::endl;
                d = TinyAD::newton_direction(g, H_proj, solver, 0.);
                dec = TinyAD::newton_decrement(d, g);

                _solve_residual = (H_proj * d + g).squaredNorm();
                
              


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
                identity_weight = std::max(identity_weight / 2., identity_min);
                if (identity_weight == identity_min)
                {
                  identity_vanished = true;
                }
                else 
                {
                  identity_vanished = false;
                }
              

                _solve_residual = (H * d + g).squaredNorm();
              }
              
            }
            catch(const std::exception& e)
            {
              try{
               std::cout << "*** diagonally regularized hessian with identity weight: " << identity_weight << " not PSD. failed to produce a descent direction.  Falling back on SVD projected hessian + increasing regularization weight. ***" << std::endl;

              this->eval_func_and_proj_hess_to_psd_local(x);
              f = this->get_fval_at_x(); // maybe don't need these two lines 
              g = this->get_grad_at_x();
              H_proj = this->get_hessian_at_x();
              d = TinyAD::newton_direction(g, H_proj, solver);
              dec = TinyAD::newton_decrement(d, g);

              _solve_residual = (H_proj * d + g).squaredNorm();
              if ( !useProjHessian )
                identity_weight = identity_weight * 10.;
              else 
                assert(false); // this should never happen
              }
              catch(const std::exception& e)
              {
                std::cout << x << std::endl;
                std::cout << H_proj << std::endl;
              }
            }

            std::cout << "Solve Residual " << std::scientific << _solve_residual << " | gradient norm " << _rhs_norm << " | ratio " << _solve_residual / _rhs_norm << std::endl;

          */



            auto t2 = std::chrono::high_resolution_clock::now();
            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            // std::cout << ms_int.count() << "ms\n";


            // _cur_x = TinyAD::line_search(x, d, f, g, [&] (Eigen::VectorXd& point) -> double { return this->eval_func_at(point); }, 1., .8, 512, 1e-3);
            _newton_dir = d;
            _dec = dec;
            _max_gradient_norm = g.cwiseAbs().maxCoeff(); 
            _prev_step_progress = f - this->eval_func_at(_cur_x);
            _prev_step_time = ms_int.count();
            // _line_search_step_size = ; // can get this from cur_x, x and d


            // also log the gradient norm 

            // _dec = f - this->eval_func_at(_cur_x);  // use this to track the true decrement 
            
            // 
            std::cout << "prev obj " << f << " | _prev_step_progress: " << _prev_step_progress << " | newton dec: " << _dec << "  | max_gradient_norm: " << _max_gradient_norm << std::endl;

          std::cout << "finshed newton step" << " | solve took " << _prev_step_time << " ms " << std::endl;


            return _cur_x;
    };