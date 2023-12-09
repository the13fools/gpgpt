#ifndef SOLVE_CURL_PROJ_H
#define SOLVE_CURL_PROJ_H



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


class Solve_Curl_Proj 
public:
    Solve_Curl_Proj() {

    
    }

    ~Solve_Curl_Proj(){
      delete _opt;
    }

    void initSimulation(AppState appState_input)
    {



      // move this inside mint2d
      appState->solveStatus = "init curl projection solver";

      std::cout << "**** setup curl projection tinyAD optimization ****" << std::endl;

      func = TinyAD::scalar_function<DOFS_PER_ELEMENT>(TinyAD::range(appState->F.rows()));

      /////////////////////////////
      /// 
      ///  setup the type of curl operator
      ///
      /// TODO make this a little more generic later
      /// 
      /////////////////////////////
      // OptZoo::addConstTestTerm(func, *appState);
      
      // OptZoo<DOFS_PER_ELEMENT>::addPinnedBoundaryTerm(func, *appState);

      // OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L2_Term(func, *appState);
      // OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L2x2_Term(func, *appState);
      // OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L4_Term(func, *appState);

      OptZoo<DOFS_PER_ELEMENT>::addCurlTerm(func, *appState);

      // Update params specific to this solve here

      appState->config->w_attenuate = 1.;
      appState->config->w_smooth = 1e5;
      appState->config->w_curl = 1e1;
      


  // This is some magic to get the tinyAD function to work with the ADFuncRunner interface.
  // Basically this file needs to store the actual instanced _opt, but then 
  // it also needs to be registered on the MintHook2D (which this inherits from)
  // In order to run the optimization.  
  //
  // 
  // This little bit of boiler plate goes a long way.  
      _opt = new ADFunc_TinyAD_Instance<DOFS_PER_ELEMENT>();
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
        // appState->deltas.row(i) = Eigen::VectorXd::Zero(4);
        x.segment<2>(nvars*i) = appState->frames.row(i);
        // x.segment<4>(nvars*i+2) = appState->deltas.row(i);
        
      }
      _opt->_cur_x = x;

      appState->config->w_smooth_vector = 0;

    }







protected:
  // Read mesh and compute Tutte embedding




    ADFunc_TinyAD_Instance<DOFS_PER_ELEMENT>* _opt;
    decltype(TinyAD::scalar_function<DOFS_PER_ELEMENT>(TinyAD::range(1))) func;



    
};





#endif // SOLVE_L2_NETWON_RANK1_H

