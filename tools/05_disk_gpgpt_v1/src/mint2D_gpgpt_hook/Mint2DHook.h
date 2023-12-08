#ifndef MINT2DHOOK_H
#define MINT2DHOOK_H

#include "AppState.h" // Include AppState definition
#include <Eigen/Dense>
#include <Eigen/Core>
#include <string>

#include "PhysicsHook.h"
#include "ADWrapper/ADFuncRunner.h"
#include "ADWrapper/ADFunc_TinyAD_Instance.h"

#include "Surface.h"
#include "FileParser.h"

#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/LinearSolver.hh>


enum class DOFType {
    primals, moments, deltas, gammas, Element_COUNT
};

class Mint2DHook : public virtual PhysicsHook
{
public:

    Mint2DHook(AppState* state) : PhysicsHook() {
        appState = state;
        outputData = appState->os;
        opt = new ADFunc_TinyAD_Instance<6>();
    //   current_element = Field_View::vec_norms;
    }
    virtual ~Mint2DHook(){
        delete appState;
        delete outputData;
        delete opt;
    }
    
    virtual void drawGUI();
    virtual void updateRenderGeometry();
    virtual void renderRenderGeometry();
    virtual void initSimulation();
    virtual bool simulateOneStep();
    virtual void pause();
    virtual void updateAppStateFromOptState(){ return; };

    void resetAppState();
    void initializeLogFolder(); 
    void initializeOtherParameters(); 
    void initBoundaryConditions();
    void initCurlOperators();


// these are moved to appState
    // Eigen::MatrixXd V; // Vertex positions
    // Eigen::MatrixXi F; // Face indices
    // Eigen::MatrixXd frames_orig;
    // Eigen::MatrixXd frames;

    // current state stored in opt
    // Eigen::VectorXd x;

    TinyAD::LinearSolver<double> solver; // make this changable 
    // decltype(TinyAD::scalar_function<6>(TinyAD::range(1))) func;
    Eigen::VectorXi bound_face_idx; // the faces on the boundary, for now let tinyAD do the boundary enforcement 


    AppState* appState; // Pointer to AppState instance

    OutputState* outputData; // threadsafe copy of the output fields.

    // AppState* renderState; // We clone visualization data to a different object, and update it in a threadsafe manner.  
    
    // This is an interface to take optimization steps with respect to a set objective.
    // there is some magic going on under the hood mostly involving the tinyAD library at the moment.
    // this wrapper is templated and can instead inject your own autodiff framework or solver instead.
    ADFuncRunner* opt;

    std::unique_ptr<FileParser> fileParser;

private:
    void updateOptimizationParameters();
    void checkAndUpdateConvergence(double decrement, double energy);
    void updateVisualizationData();
    void finalizeIteration();




    // Other private member variables and functions as needed
};

#endif // MINT2DHOOK_H
