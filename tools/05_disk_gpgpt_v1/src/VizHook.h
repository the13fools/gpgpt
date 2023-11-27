#include "PhysicsHook.h"
#include "Mint2DHook.h"

#include "Surface.h"


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>

#include "ADWrapper/ADFuncRunner.h"
#include "ADWrapper/ADFunc_TinyAD_Instance.h"



#include <igl/writeDMAT.h>

#include <sys/stat.h>

#include <igl/map_vertices_to_circle.h>



#include <chrono>

// #include <fstream>
#include <sys/stat.h>


#include "OptZoo.h"



class VizHook : public Mint2DHook
{
public:
    VizHook() : Mint2DHook(new AppState()) {
      appState->current_element = Field_View::vec_norms;
    }


    virtual void drawGUI()
    {
      Mint2DHook::drawGUI();

    }

    ~VizHook(){
      delete _opt;
    }




    virtual void initSimulation()
    {

      // cur_mesh_name = "circle_subdiv";

      // cur_mesh_name = "circle";
      appState->meshName = "circle_irreg";
      // appState->meshName = "circle_1000";



      // cur_mesh_name = "circle_irreg_20000";


      // cur_mesh_name = "circle_1000";


      // Call Parent initialization to load mesh and initialize data structures
      // Add file parsing logic here.
      Mint2DHook::initSimulation();


      std::cout << "**** setup tinyAD optimization ****" << std::endl;

      // In this example we are setting up an optimization with a vector and a delta per face.  

//       // Set up function with 2D vertex positions as variables.
//       func = TinyAD::scalar_function<6>(TinyAD::range(F.rows()));

    func = TinyAD::scalar_function<6>(TinyAD::range(appState->F.rows()));
    // auto func = _opt._func;

    // setup tinyad func


    /////////////////////////////
    /// 
    /////////////////////////////
    OptZoo::addConstTestTerm(func, *appState);
    OptZoo::addPinnedBoundaryTerm(func, *appState);

    // OptZoo::addPinnedBoundaryTerm(func, *appState);
    OptZoo::addSmoothnessTerm(func, *appState);
    


// This is some magic to get the tinyAD function to work with the ADFuncRunner interface.
    _opt = new ADFunc_TinyAD_Instance<6>();
    _opt->set_tinyad_objective_func(&func);
    // _opt->_func = &func;
    opt = static_cast<ADFuncRunner*>(_opt);

    std::cout << "DOFS in opt" << opt->_cur_x.rows() << std::endl;
    // std::cout << "nvars in opt" << _opt->_func->n_vars << std::endl; // get_num_vars
    std::cout << "nvars in opt" << this->opt->get_num_vars() << std::endl; 




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

    virtual void updateAppStateFromOptState()
    {
      int nelem = appState->F.rows();
      int nvars = 6; // opt->get_num_vars();

      Eigen::VectorXd x = opt->get_current_x();
      for(int i = 0; i < nelem; i++)
      {
        appState->frames.row(i) = x.segment<2>(nvars*i);
        appState->deltas.row(i) = x.segment<4>(nvars*i+2);

        
      }

        // std::cout << x.segment<2>(nvars*i) << std::endl;
        // std::cout << x.segment<2>(nvars*i+4) << std::endl;

        // std::cout << appState->frames.rows() << " " << appState->deltas.rows() << std::endl;
        // std::cout << appState->frames.cols() << " " << appState->deltas.cols() << std::endl;


      // std::cout << "updateAppStateFromOptState" << std::endl;
      //   func.x_to_data(opt->get_current_x(), [&] (int f_idx, const Eigen::VectorXd& v) {
      //           appState->frames.row(f_idx) = v.head<2>();
      //           appState->deltas.row(f_idx) = v.tail<2>();
      //           // if (bound_face_idx(f_idx) == 1)
      //           // {
      //           //   frames.row(f_idx) = frames_orig.row(f_idx);
      //           // }
      //           });
    }




protected:
  // Read mesh and compute Tutte embedding



  // Eigen::MatrixXd metadata;

     ADFunc_TinyAD_Instance<6>* _opt;
     decltype(TinyAD::scalar_function<6>(TinyAD::range(1))) func;


      std::vector<Eigen::Matrix2d> rots;// 
      std::vector<Eigen::Matrix4d> rstars;
      std::vector<Eigen::Vector4d> e_projs;
      std::vector<Eigen::Vector2d> e_projs_primal;

      Eigen::MatrixXd e_projs2;



  


  // Eigen::MatrixXd renderFrames;
  // Eigen::MatrixXd renderDeltas;





  
  std::vector<Eigen::Matrix2d> rest_shapes;
// %%% 2 + 4 + 4







  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg_solver;


    
};


// TODO initialize app state.  

//       // Assemble inital x vector from P matrix.
//       // x_from_data(...) takes a lambda function that maps
//       // each variable handle (vertex index) to its initial 2D value (Eigen::Vector2d).
//         x = func.x_from_data([&] (int f_idx) {
//           Eigen::VectorXd ret;
//           ret = Eigen::VectorXd::Zero(6); // resize(10);
//           ret.head(2) = Eigen::VectorXd::Random(2) * 0.0001;
//           // ret << frames.row(f_idx), deltas.row(f_idx);
//           return ret;
//           });