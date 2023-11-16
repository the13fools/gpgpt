#ifndef MINT2DHOOK_H
#define MINT2DHOOK_H

#include "AppState.h" // Include AppState definition
#include <Eigen/Dense>
#include <string>

#include "PhysicsHook.h"

#include "Surface.h"

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/LinearSolver.hh>

class Mint2DHook : public virtual PhysicsHook
{
public:

    Mint2DHook(AppState* state) : PhysicsHook() {
        appState = state;
    //   current_element = Field_View::vec_norms;
    }
    virtual ~Mint2DHook(){
        delete appState;
        delete renderState;
        delete opt;
    }
    
    virtual void drawGUI();
    virtual void updateRenderGeometry();
    virtual void renderRenderGeometry();
    virtual void initSimulation();
    virtual bool simulateOneStep();

    // void reset();
    void initializeLogFolder(); 
    void initializeOtherParameters(); 
    void initBoundaryConditions();


    std::string cur_mesh_name;
    Surface cur_surf;
    Eigen::MatrixXd V; // Vertex positions
    Eigen::MatrixXi F; // Face indices
    Eigen::MatrixXd frames_orig;
    Eigen::MatrixXd frames;
    Eigen::VectorXd x;

    TinyAD::LinearSolver<double> solver; // make this changable 
    decltype(TinyAD::scalar_function<6>(TinyAD::range(1))) func;
    Eigen::VectorXi bound_face_idx; // the faces on the boundary, for now let tinyAD do the boundary enforcement 


    AppState* appState; // Pointer to AppState instance
    AppState* renderState; // We clone visualization data to a different object, and update it in a threadsafe manner.  
    
    // This is an interface to take optimization steps with respect to a set objective.
    // there is some magic going on under the hood mostly involving the tinyAD library at the moment.
    // this wrapper is templated and can instead inject your own autodiff framework or solver instead.
    ADFuncRunner* opt; 

private:
    void updateOptimizationParameters();
    void checkAndUpdateConvergence(double decrement, double energy);
    void updateVisualizationData();
    void finalizeIteration();




    // Other private member variables and functions as needed
};

#endif // MINT2DHOOK_H
