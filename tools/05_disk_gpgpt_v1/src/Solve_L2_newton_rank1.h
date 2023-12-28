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

#include <igl/on_boundary.h>

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
      // delete _opt;
    }

    virtual void drawGUI()
    {
      Mint2DHook::drawGUI();

    }

    virtual void initSimulation()
    {


        // appState->meshName = "circle_1000";
      // appState->meshName = "circle_subdiv";
      // appState->meshName = "circle";
      // appState->meshName = "circle_irreg";
      // appState->meshName = "circle_irreg_20000";
      appState->meshName = "disk_irreg";
      


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

      OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L2_Term(func, *appState);
      // OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L2x2_Term(func, *appState);
      OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L4_Term(func, *appState);

      OptZoo<DOFS_PER_ELEMENT>::addCurlTerm_L2(func, *appState);
      OptZoo<DOFS_PER_ELEMENT>::addCurlTerm_L4(func, *appState);


      // Update params specific to this solve here
      


  // This is some magic to get the tinyAD function to work with the ADFuncRunner interface.
  // Basically this file needs to store the actual instanced _opt, but then 
  // it also needs to be registered on the MintHook2D (which this inherits from)
  // In order to run the optimization.  
  //
  // 
  // This little bit of boiler plate goes a long way.  
      _opt = std::make_unique< ADFunc_TinyAD_Instance<DOFS_PER_ELEMENT> >();
      _opt->set_tinyad_objective_func(&func);
      if(appState->shouldReload)
      {
        updateOptStateFromAppState();
        appState->shouldReload = false;
      }
      else
      {
        init_opt_state();
      }
      // _opt->_func = &func;
      opt = static_cast<ADFuncRunner*>(_opt.get());

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
        // appState->frames.row(i) = Eigen::VectorXd::Random(2) * 1e-1;

        // -4 r Cos[(\[Theta] - Pi)/2], r Sin[(\[Theta] - Pi)/2]

        Eigen::RowVector3d centroid = (appState->V.row(appState->F(i, 0)) +
                                            appState->V.row(appState->F(i, 1)) +
                                            appState->V.row(appState->F(i, 2))) / 3.0;

        double r = centroid.norm() + 1e-10;
        double theta = atan2(centroid(1),centroid(0)); // acos(vec(1)) * .5;

        // Eigen::Vector2d frame_cur = 1./r * Eigen::Vector2d(centroid.y(), -centroid.x()).normalized();
        Eigen::Vector2d frame_cur = r * Eigen::Vector2d( -4 * ( cos( (theta - 3.1415)  / 2. ) ), r * sin( (theta - 3.1415) / 2. ) ).normalized();

        appState->frames.row(i) = frame_cur;

  



        // appState->frames.row(i)

      // appState->problemFileTag = "_+0_5_";
      // frame(0) = cos(theta * .5);
      // frame(1) = sin(theta * .5);
 
        // appState->deltas.row(i) = Eigen::VectorXd::Zero(4);
        if (appState->bound_face_idx(i) == 1) {
          appState->frames.row(i) = appState->frames_orig.row(i);
        }
        x.segment<2>(nvars*i) = appState->frames.row(i);
        // x.segment<4>(nvars*i+2) = appState->deltas.row(i);
        
      }
      _opt->_cur_x = x;

      appState->config->w_smooth_vector = 0;

    }

    virtual void initBoundaryConditions() {
        // Assuming boundary faces are identified in AppState
        Eigen::MatrixXi K;

        // Eigen::MatrixXi bound_face_idx = appState->bound_face_idx;

        Eigen::VectorXi boundaryFaces;
        igl::on_boundary(appState->F,boundaryFaces, K);

        appState->bound_face_idx = boundaryFaces;

        appState->frames.resize(appState->F.rows(), DOFS_PER_ELEMENT);

        // Initialize boundary conditions
        for (int i = 0; i < boundaryFaces.size(); ++i) {
            if (boundaryFaces(i) == 1) { // If face is on the boundary
                Eigen::RowVector3d centroid = (appState->V.row(appState->F(i, 0)) +
                                            appState->V.row(appState->F(i, 1)) +
                                            appState->V.row(appState->F(i, 2))) / 3.0;

                if (centroid.norm() < 0.45) { // Custom condition for boundary faces
                    boundaryFaces(i) = -1; // Mark for special handling or exclusion
                } else {
                    // Set frame orientation based on the centroid

                    // 
                    Eigen::Vector2d vec = Eigen::Vector2d(centroid.x(), centroid.y()).normalized();

                    Eigen::VectorXd frame = Eigen::VectorXd::Zero(2);
                    double theta = atan2(vec(1),vec(0)); // acos(vec(1)) * .5;


                    // -4 r Cos[(\[Theta] - Pi)/2], r Sin[(\[Theta] - Pi)/2]


                    // appState->problemFileTag = "_radial_";
                    // frame(0) = cos(theta);
                    // frame(1) = sin(theta);

                    appState->problemFileTag = "_+0_5_";
                    frame(0) = cos(theta * .5);
                    frame(1) = sin(theta * .5);


                    // appState->problemFileTag = "_-0_5_";
                    // frame(0) = -4 * cos((theta-3.1415) * .5 );
                    // frame(1) = sin((theta-3.1415) * .5 );
                    // frame(2) = -vec(1);
                    // frame(3) = vec(0);
                    appState->frames.row(i) = frame;
                }
            }
        }

        appState->frames_orig = appState->frames;

        // 

        // 

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
      appState->config->w_smooth = 1e0;
      appState->config->w_bound = 1e5;
      appState->config->w_curl = 1e-3;
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
      if ( appState->keepSolving == false && appState->config->w_attenuate > 1e-14)
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




    std::unique_ptr<ADFunc_TinyAD_Instance<DOFS_PER_ELEMENT> > _opt;
    decltype(TinyAD::scalar_function<DOFS_PER_ELEMENT>(TinyAD::range(1))) func;



    
};





#endif // SOLVE_L2_NETWON_RANK1_H

