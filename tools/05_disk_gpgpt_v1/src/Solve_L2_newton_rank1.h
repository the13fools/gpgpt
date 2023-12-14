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

#define DOFS_PER_ELEMENT 2


class Solve_L2_newton_rank1 : public Mint2DHook
{
public:
    Solve_L2_newton_rank1() : Mint2DHook(new AppState()) {
      appState->current_element = Field_View::vec_norms;
      appState->solveType = "L2_newton_rank1";
      appState->solveDescription = "L2_newton_rank1";




      appState->primals_layout = {0, 2};
      appState->moments_layout = {0, 0};
      // appState->deltas_layout = {2, 4};
      appState->deltas_layout = {2, 0};

      assert(DOFS_PER_ELEMENT == (appState->primals_layout.size + appState->moments_layout.size + appState->deltas_layout.size));

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

      // appState->meshName = "circle_1000";
      // appState->meshName = "circle_subdiv";
      appState->meshName = "circle";
      // appState->meshName = "circle_irreg";
      // appState->meshName = "circle_irreg_20000";
      


      // Call Parent initialization to load mesh and initialize data structures
      // Add file parsing logic here.
      Mint2DHook::initSimulation();

      // move this inside mint2d
      appState->solveStatus = "init L2 newton rank 1";

      std::cout << "**** setup tinyAD optimization ****" << std::endl;

      func = TinyAD::scalar_function<DOFS_PER_ELEMENT>(TinyAD::range(appState->F.rows()));

      /////////////////////////////
      /// Add terms to the objective function here.  
      /// Feel free to implement your own terms :). 
      ///
      /// There's a lot of opportunity for both human and machine learning 
      /// by exploring modifications to these operators.  
      ///
      /// convenient auto-diff with sparse hessians is quite new, there's a lot to explore! 
      /////////////////////////////
      // OptZoo::addConstTestTerm(func, *appState);
      
      OptZoo<DOFS_PER_ELEMENT>::addPinnedBoundaryTerm(func, *appState);

      // OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L2_Term(func, *appState);
      // OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L2x2_Term(func, *appState);
      OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L4_Term(func, *appState);

      OptZoo<DOFS_PER_ELEMENT>::addCurlTerm(func, *appState);

      // Update params specific to this solve here
      


  // This is some magic to get the tinyAD function to work with the ADFuncRunner interface.
  // Basically this file needs to store the actual instanced _opt, but then 
  // it also needs to be registered on the MintHook2D (which this inherits from)
  // In order to run the optimization.  
  //
  // 
  // This little bit of boiler plate goes a long way.  
      _opt = new ADFunc_TinyAD_Instance<DOFS_PER_ELEMENT>();
      _opt->set_tinyad_objective_func(&func);
      if(appState->shouldReload)
      {
        updateOptStateFromAppState();
      }
      else
      {
        init_opt_state();
      }
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
        // appState->deltas.row(i) = Eigen::VectorXd::Zero(4);
        x.segment<2>(nvars*i) = appState->frames.row(i);
        // x.segment<4>(nvars*i+2) = appState->deltas.row(i);
        
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

    virtual void initConfigValues()
    {
      appState->config->w_attenuate = 1.;
      appState->config->w_smooth = 1e5;
      appState->config->w_bound = 1e8;
      appState->config->w_curl = 1e1;
    }

// This is called after each step.  
    virtual void updateAppStateFromOptState()
    {
      int nelem = appState->F.rows();
      int nvars = DOFS_PER_ELEMENT; // opt->get_num_vars();

      // Try to add post projection curl operator here.  


      Eigen::VectorXd x = opt->get_current_x();
      for(int i = 0; i < nelem; i++)
      {
        appState->frames.row(i) = x.segment<2>(nvars*i);
        // appState->deltas.row(i) = x.segment<4>(nvars*i+2);

        
      }


      

      // Make this more generic like first write a set of configs to the outdirectory and make this advance to the next one when keepSolving is false.
      if ( appState->keepSolving == false && appState->config->w_attenuate > 1e-12)
      {
        appState->config->w_attenuate = appState->config->w_attenuate / 2.;
        appState->keepSolving = true;  
        appState->outerLoopIteration += 1;
        std::cout << "~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ attenuate set to: " << appState->config->w_attenuate << " ~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~" << std::endl;

      }



    }

    // this is called after reloading data in order to be able to take forward steps.
    virtual void updateOptStateFromAppState()
    {
      Eigen::VectorXd x = _opt->_cur_x;
      int nelem = appState->F.rows();
      int nvars = DOFS_PER_ELEMENT; // opt->get_num_vars();

      // Eigen::VectorXd x = opt->get_current_x();
      for(int i = 0; i < nelem; i++)
      {
        x.segment<2>(nvars*i) = appState->frames.row(i);
      }
      _opt->_cur_x = x;

      // appState->config->w_smooth_vector = 0;

    }




protected:
  // Read mesh and compute Tutte embedding




    ADFunc_TinyAD_Instance<DOFS_PER_ELEMENT>* _opt;
    decltype(TinyAD::scalar_function<DOFS_PER_ELEMENT>(TinyAD::range(1))) func;



    
};





#endif // SOLVE_L2_NETWON_RANK1_H

