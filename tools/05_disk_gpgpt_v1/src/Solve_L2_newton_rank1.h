#ifndef SOLVE_L2_NETWON_RANK1_H
#define SOLVE_L2_NETWON_RANK1_H



#include "PhysicsHook.h"
#include "Mint2DHook.h"

#include "Surface.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <TinyAD/ScalarFunction.hh>

#include "ADWrapper/ADFuncRunner.h"
#include "ADWrapper/ADFunc_TinyAD_Instance.h"



#include <igl/writeDMAT.h>

#include "OptZoo.h"



class Solve_L2_newton_rank1 : public Mint2DHook
{
public:
    Solve_L2_newton_rank1() : Mint2DHook(new AppState()) {
      appState->current_element = Field_View::vec_norms;
      appState->solveType = "L2_newton_rank1";
      appState->solveDescription = "L2_newton_rank1";


      appState->primals_layout = {0, 2};
      appState->moments_layout = {0, 0};
      appState->deltas_layout = {2, 4};
    }

    ~Solve_L2_newton_rank1(){
      delete _opt;
    }

    virtual void drawGUI()
    {
      Mint2DHook::drawGUI();

    }

    virtual void initSimulation()
    {

      appState->meshName = "circle_subdiv";

      // appState->meshName = "circle";
      // appState->meshName = "circle_irreg";
      // appState->meshName = "circle_1000";



      // appState->meshName = "circle_irreg_20000";


      // appState->meshName = "circle_1000";
      


      // Call Parent initialization to load mesh and initialize data structures
      // Add file parsing logic here.
      Mint2DHook::initSimulation();

      appState->solveStatus = "init L2 newton rank 1";


      std::cout << "**** setup tinyAD optimization ****" << std::endl;

      // In this example we are setting up an optimization with a vector and a delta per face.  

//       // Set up function with 2D vertex positions as variables.
//       func = TinyAD::scalar_function<6>(TinyAD::range(F.rows()));

    func = TinyAD::scalar_function<6>(TinyAD::range(appState->F.rows()));
    // auto func = _opt._func;

    // setup tinyad func


    /////////////////////////////
    /// Add terms to the objective function here.  
    /// Feel free to implement your own terms :). 
    ///
    /// There's a lot of opportunity for both human and machine learning 
    /// by exploring modifications to these operators.  
    ///
    /// Auto-diff with sparse hessians is very new, there's a lot to explore! 
    /////////////////////////////
    // OptZoo::addConstTestTerm(func, *appState);
    OptZoo::addPinnedBoundaryTerm(func, *appState);

    // OptZoo::addPinnedBoundaryTerm(func, *appState);
    OptZoo::addSmoothnessTerm(func, *appState);
    OptZoo::addCurlTerm(func, *appState);


    appState->config->w_attenuate = 1.;
    appState->config->w_smooth = 1e5;
    appState->config->w_curl = 1e1;
    


// This is some magic to get the tinyAD function to work with the ADFuncRunner interface.
    _opt = new ADFunc_TinyAD_Instance<6>();
    _opt->set_tinyad_objective_func(&func);
    init_opt_state();
    // _opt->_func = &func;
    opt = static_cast<ADFuncRunner*>(_opt);

    std::cout << "DOFS in opt" << opt->_cur_x.rows() << std::endl;
    // std::cout << "nvars in opt" << _opt->_func->n_vars << std::endl; // get_num_vars
    std::cout << "nvars in opt" << this->opt->get_num_vars() << std::endl; 




    }

// Init state and clone it to the appstate in order to make the visualization accurate.  
    void init_opt_state()
    {
      Eigen::VectorXd x = _opt->_cur_x;
      int nelem = appState->F.rows();
      int nvars = DOFS_PER_ELEMENT; // opt->get_num_vars();

      // Eigen::VectorXd x = opt->get_current_x();
      for(int i = 0; i < nelem; i++)
      {
        appState->frames.row(i) = Eigen::VectorXd::Random(2) * 1e-1;
        appState->deltas.row(i) = Eigen::VectorXd::Zero(4);
        x.segment<2>(nvars*i) = appState->frames.row(i);
        x.segment<4>(nvars*i+2) = appState->deltas.row(i);
        
      }
      _opt->_cur_x = x;

      appState->config->w_smooth_vector = 0;

    }


    virtual void updateRenderGeometry()
    {
      Mint2DHook::updateRenderGeometry();

    }




    virtual void renderRenderGeometry()
    {
      Mint2DHook::renderRenderGeometry();
    }


    virtual bool simulateOneStep()
    {
      return Mint2DHook::simulateOneStep();

    }

// This is called after each step.  
    virtual void updateAppStateFromOptState()
    {
      int nelem = appState->F.rows();
      int nvars = DOFS_PER_ELEMENT; // opt->get_num_vars();

      Eigen::VectorXd x = opt->get_current_x();
      for(int i = 0; i < nelem; i++)
      {
        appState->frames.row(i) = x.segment<2>(nvars*i);
        appState->deltas.row(i) = x.segment<4>(nvars*i+2);

        
      }


      

            // Make this more generic like first write a set of configs to the outdirectory and make this advance to the next one when keepSolving is false.
      if ( appState->keepSolving == false && appState->config->w_attenuate > 1e-12)
      {
        appState->config->w_attenuate = appState->config->w_attenuate / 10.;
        appState->keepSolving = true;  
        std::cout << "~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ attenuate set to: " << appState->config->w_attenuate << " ~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~" << std::endl;

      }


      // // Make this more generic like first write a set of configs to the outdirectory and make this advance to the next one when keepSolving is false.
      // if ( appState->keepSolving == false && appState->config->w_smooth_vector > 0)
      // {
      //   appState->keepSolving = true;  
      //   appState->config->w_smooth_vector = 0;
      //   opt->useProjHessian = true; // reset to use PSD hessian because optimization problem changed.

      //   std::cout << "~~~~~~switch off primal smoothness term used to initialize the opt~~~~~" << std::endl;

      // }
      // std::cout << "appState->config->w_smooth_vector " << appState->config->w_smooth_vector << std::endl;

    }




protected:
  // Read mesh and compute Tutte embedding



  // Eigen::MatrixXd metadata;
    int DOFS_PER_ELEMENT = 6;
    ADFunc_TinyAD_Instance<6>* _opt;
    decltype(TinyAD::scalar_function<6>(TinyAD::range(1))) func;




  
  std::vector<Eigen::Matrix2d> rest_shapes;
// %%% 2 + 4 + 4


  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;


    
};





#endif // SOLVE_L2_NETWON_RANK1_H