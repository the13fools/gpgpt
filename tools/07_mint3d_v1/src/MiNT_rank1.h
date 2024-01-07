#ifndef MINT_RANK1_H
#define MINT_RANK1_H



#include "PhysicsHook.h"
#include "Mint3DHook.h"

#include "Surface.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <TinyAD/ScalarFunction.hh>

#include "ADWrapper/ADFuncRunner.h"
#include "ADWrapper/ADFunc_TinyAD_Instance.h"

// #include <igl/on_boundary.h>

// #include <igl/writeDMAT.h>

#include "OptZoo.h"

#define DOFS_PER_ELEMENT 9


class MiNT_rank1 : public Mint3DHook
{
public:
    MiNT_rank1() : Mint3DHook(new AppState()) {
      appState->current_element = Field_View::vec_norms;
      appState->solveType = "MiNT_rank1";
      appState->solveDescription = "MiNT_rank1";




      appState->primals_layout = {0, DOFS_PER_ELEMENT, DOFS_PER_ELEMENT / 3};
      appState->moments_layout = {0, 0};
      // appState->deltas_layout = {2, 4};
      appState->deltas_layout = {3, 0};

      assert(DOFS_PER_ELEMENT == (appState->primals_layout.size + appState->moments_layout.size + appState->deltas_layout.size));

      appState->useBoundaryFrames = true;

    }

    ~MiNT_rank1(){
      // delete _opt;
    }

    virtual void initConfigValues()
    {
      appState->config->w_attenuate = 1.;
      appState->config->w_smooth = 1e0;
      appState->config->w_bound = 1e1;
      appState->config->w_curl = 7e-10; // this is off for the first outer iter 
    }

    virtual void drawGUI()
    {
      Mint3DHook::drawGUI();

    }

    virtual void initSimulation()
    {


        // appState->meshName = "circle_1000";
      // appState->meshName = "circle_subdiv";
      // appState->meshName = "circle";
      // appState->meshName = "circle_irreg";
      // appState->meshName = "circle_irreg_20000";

      // appState->meshName = "disk_v210";
      // appState->meshName = "disk_v623";
      // appState->meshName = "disk_v1000";
      // appState->meshName = "disk_3480_tets";
  // appState->meshName = "cylinder_400";
  // appState->meshName = "cylinder3k";

  // appState->meshName = "sphere_r0.17";
    // appState->meshName = "sphere_r0.14";
    // appState->meshName = "sphere_r0.10";
        // appState->meshName = "sphere_r0.05";

// appState->meshName = "triangle_notwist_400";
appState->meshName = "tetrahedron_200";
      // tetrahedron_100
// parallelogram_exact

// appState->meshName = "parallelogram_exact";
// appState->meshName = "tetrahedron";

      


      // Call Parent initialization to load mesh and initialize data structures
      // Add file parsing logic here.
      Mint3DHook::initSimulation();

      // move this inside mint2d
      appState->solveStatus = "init MiNT rank1";

      std::cout << "**** setup tinyAD optimization ****" << std::endl;

      func = TinyAD::scalar_function<DOFS_PER_ELEMENT>(TinyAD::range(appState->nelem));


      // func = TinyAD::scalar_function<DOFS_PER_ELEMENT>(TinyAD::range(appState->F.rows()));

      /////////////////////////////
      /// Add terms to the objective function here.  
      /// Feel free to implement your own terms :). 
      ///
      /// There's a lot of opportunity for both human and machine learning 
      /// by exploring modifications to these operators.  
      ///
      /// convenient auto-diff with sparse hessians is quite new, there's a lot to explore! 
      /////////////////////////////
      // OptZoo<DOFS_PER_ELEMENT>::addConstTestTerm(func, *appState);

      OptZoo<DOFS_PER_ELEMENT>::addWeakOrthogonalityTerm(func, *appState);
      OptZoo<DOFS_PER_ELEMENT>::addUnitNormTerm(func, *appState);
      
      OptZoo<DOFS_PER_ELEMENT>::addNormalBoundaryTerm(func, *appState);
      
      
      // OptZoo<DOFS_PER_ELEMENT>::addPinnedBoundaryTerm(func, *appState);

      OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L2_Term(func, *appState);
      // OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L2x2_Term(func, *appState);
      OptZoo<DOFS_PER_ELEMENT>::addSmoothness_L4_Term(func, *appState);

      appState->curl_orders = {2,4,6};
      OptZoo<DOFS_PER_ELEMENT>::addCurlTerms(func, *appState);

      
      // OptZoo<DOFS_PER_ELEMENT>::addCurlTerm_L2(func, *appState);
      // OptZoo<DOFS_PER_ELEMENT>::addCurlTerm_L4(func, *appState);


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
      int nelem = appState->nelem;
      int nvars = DOFS_PER_ELEMENT; // opt->get_num_vars();
      int ntets = appState->T.rows();

      // Eigen::VectorXd x = opt->get_current_x();
      for(int i = 0; i < ntets; i++)
      {
        // appState->frames.row(i) = Eigen::VectorXd::Random(DOFS_PER_ELEMENT) * 1e-1;

        Eigen::Vector3d a = Eigen::VectorXd::Random(3); 
        Eigen::Vector3d b = Eigen::VectorXd::Random(3); 
        Eigen::Vector3d c = a.cross(b);
        b = c.cross(a);
        a.normalize();
        b.normalize();
        c.normalize();
        appState->frames.row(i).segment<3>(0) = a;
        if (DOFS_PER_ELEMENT >= 6)
        {
          appState->frames.row(i).segment<3>(3) = b;
        }
        if (DOFS_PER_ELEMENT >= 9)
        {
          appState->frames.row(i).segment<3>(6) = c;
        }


        // for (int j = 0; j < DOFS_PER_ELEMENT / 3; ++j)
        // {
        //   appState->frames.row(i).segment<3>(j*3) = 0.;
        // }
       
       
       
        appState->frames.row(i) = Eigen::VectorXd::Random(DOFS_PER_ELEMENT) * 1e-4;
        // appState->frames.row(i) = Eigen::VectorXd::Ones(DOFS_PER_ELEMENT);

// add back pinned constraints later.  
        // if (appState->bound_face_idx(i) == 1) {
        //   appState->frames.row(i) = appState->frames_orig.row(i);
        // }

        x.segment<DOFS_PER_ELEMENT>(nvars*i) = appState->frames.row(i);
      }

      for (int i = 0; i < nelem - ntets; ++i)
      {
        appState->boundary_frames.row(i) = Eigen::VectorXd::Random(DOFS_PER_ELEMENT) * 1e-1;
        x.segment<DOFS_PER_ELEMENT>(nvars*(i+ntets)) = appState->boundary_frames.row(i);
      }

      // add pin here 

      _opt->_cur_x = x;

      appState->config->w_smooth_vector = 0;

    }

    virtual void initBoundaryConditions() {
        // Assuming boundary faces are identified in AppState
        Eigen::MatrixXi K;

        appState->frames.resize(appState->T.rows(), DOFS_PER_ELEMENT);
        appState->boundary_frames.resize(appState->cur_tet_mesh->nBoundaryElements(), DOFS_PER_ELEMENT);


        // Eigen::MatrixXi bound_face_idx = appState->bound_face_idx;
/*
        Eigen::VectorXi boundaryFaces;
        igl::on_boundary(appState->T,boundaryFaces, K);

        appState->bound_face_idx = boundaryFaces;


        // Initialize boundary conditions
        for (int i = 0; i < boundaryFaces.size(); ++i) {
            if (boundaryFaces(i) == 1) { // If face is on the boundary
                Eigen::RowVector3d centroid = (appState->V.row(appState->T(i, 0)) +
                                            appState->V.row(appState->T(i, 1)) +
                                            appState->V.row(appState->T(i, 2)) + 
                                            appState->V.row(appState->T(i, 3))) / 4.0;

                centroid(2) = 0.;

                // if (centroid.norm() < 0.45) { // Custom condition for boundary faces
                //     boundaryFaces(i) = -1; // Mark for special handling or exclusion
                // if (centroid.norm() < .3) { // Custom condition for boundary faces
                // if (centroid.norm() < 4) { // disk_3480_tets
                // if (centroid.norm() < 45) { // disk_v623
                // if (centroid.norm() < 95) { // cylinder3k
                if (centroid.norm() < .2) { // sphere



                    boundaryFaces(i) = 0;
                } else {
                    // Set frame orientation based on the centroid
                    Eigen::Vector2d vec = Eigen::Vector2d(centroid.x(), centroid.y()).normalized();

                    Eigen::VectorXd frame = Eigen::VectorXd::Zero(DOFS_PER_ELEMENT);
                    double theta = atan2(vec(1),vec(0)); // acos(vec(1)) * .5;
                    // frame(0) = cos(theta);
                    // frame(1) = sin(theta);
                    // frame(2) = 0.;


                    frame(0) = -sin(theta);
                    frame(1) = cos(theta);
                    frame(2) = 0.;
                    // frame(2) = -vec(1);
                    // frame(3) = vec(0);
                    appState->frames.row(i) = frame;

                    // appState->frames.row(i) = frame *  1./(centroid.norm() + 1e-10);
                    // appState->frames.row(i) = frame * (centroid.norm() + 1e-10);

                }
            }
        }

        appState->bound_face_idx = boundaryFaces;
        appState->frames_orig = appState->frames;
 
*/

        // 

        // 

    }


    virtual void updateRenderGeometry()
    {
      Mint3DHook::updateRenderGeometry();

    }




    virtual void renderRenderGeometry()
    {
      Mint3DHook::renderRenderGeometry();
    }


    virtual bool simulateOneStep()
    {
      return Mint3DHook::simulateOneStep();

    }



// This is called after each step.  
    virtual void updateAppStateFromOptState()
    {
      int nelem = appState->nelem;
      int nvars = DOFS_PER_ELEMENT; // opt->get_num_vars();
      int ntets = appState->T.rows();

      // Try to add post projection curl operator here.  


      Eigen::VectorXd x = opt->get_current_x();
      for(int i = 0; i < ntets; i++)
      {
        appState->frames.row(i) = x.segment<DOFS_PER_ELEMENT>(nvars*i);

        // if (appState->bound_face_idx(i) == 1) {
        //   appState->frames.row(i) = appState->frames_orig.row(i);
        // }
        // appState->deltas.row(i) = x.segment<4>(nvars*i+2);
    
        
      }


      for(int i = 0; i < nelem - ntets; i++)
      {
        appState->boundary_frames.row(i) = x.segment<DOFS_PER_ELEMENT>( nvars*(i + ntets) );   
        
      }


      if ( appState->keepSolving == false && 1./appState->config->w_attenuate * appState->config->w_curl < 1e17)
      {
        appState->config->w_attenuate = appState->config->w_attenuate / 2.;
        appState->keepSolving = true;  
        appState->outerLoopIteration += 1;
        std::cout << "~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ attenuate set to: " << appState->config->w_attenuate << " ~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~" << std::endl;

      }
      else if ( appState->keepSolving == false && 1./appState->config->w_attenuate * appState->config->w_curl < 1e23)
      {
        appState->config->w_attenuate = appState->config->w_attenuate / 10.;
        appState->keepSolving = true;  
        appState->outerLoopIteration += 1;
        std::cout << "~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ attenuate set to: " << appState->config->w_attenuate << " ~~~~~ ~~~~~~ ~~~~~~ ~~~~~~ ~~~~~~" << std::endl;
      }




    }

    // this is called after reloading data in order to be able to take forward steps.
    virtual void updateOptStateFromAppState()
    {
      Eigen::VectorXd x = _opt->_cur_x;
      int nelem = appState->nelem;
      int ntets = appState->T.rows();
      int nvars = DOFS_PER_ELEMENT; // opt->get_num_vars();

      // Eigen::VectorXd x = opt->get_current_x();
      for(int i = 0; i < ntets; i++)
      {
        x.segment<DOFS_PER_ELEMENT>(nvars*i) = appState->frames.row(i);
      }

      for (int i = 0; i < nelem - ntets; ++i)
      {
        x.segment<DOFS_PER_ELEMENT>(nvars*(i+ntets)) = appState->boundary_frames.row(i);
      }

      

      _opt->_cur_x = x;

      // appState->config->w_smooth_vector = 0;

    }




protected:
  // Read mesh and compute Tutte embedding




    std::unique_ptr<ADFunc_TinyAD_Instance<DOFS_PER_ELEMENT> > _opt;
    decltype(TinyAD::scalar_function<DOFS_PER_ELEMENT>(TinyAD::range(1))) func;



    
};





#endif // MINT_RANK1_H

#undef DOFS_PER_ELEMENT

