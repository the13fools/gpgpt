#ifndef APPSTATE_H
#define APPSTATE_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <optional>
#include <unordered_map>
#include <limits>

#include "MyConfig.h"
#include "FieldView.h"

#include "Surface.h"
#include <memory>
 
using Field_View = Views::Field_View;

// Enum for identifying field view quantities
// enum Field_View {
//     vec_norms, delta_norms, vec_dirch, moment_dirch, primal_curl_residual, sym_curl_residual, gui_free, Element_COUNT
// };
// Views::

// enum Field_View : unsigned int;

// Struct to hold bounds for each field view quantity
struct FieldBounds {
    float upper = .9; // std::numeric_limits<float>::max();
    float lower = .1; //std::numeric_limits<float>::lowest();

    // std::pair<double, double> getBoundsPair() {
    //     return std::make_pair(lower, upper);
    // }
};

struct DOFMemoryLayout {
    int start = 0;
    int size = 2; 
    int rank = 1; // This is mostly used when optimizing directly over primal variables in higher rank setting.  
    // TODO implement support for this! 
};

// enum VariableType {
//     frames, moments, deltas, gammas, Element_COUNT
// };

class OutputState 
{
public: 

    std::vector<Eigen::MatrixXd> frames;

    // derived variables 
    double cur_global_objective_val;
    // Eigen::MatrixXd renderFrames;
    Eigen::VectorXd norms_vec;
    Eigen::VectorXd norms_delta;
    Eigen::VectorXd norms_moment; // TODO IMPLEMENT THIS!
    Eigen::VectorXd thetas; // TODO  IMPORTANT! 
    Eigen::VectorXd curls_primal;

    // Current curl component of the energy, this should be a constraint that 
    // is close to numerically zero at convergence.  
    Eigen::VectorXd curls_sym;
    Eigen::VectorXd smoothness_primal;

    // Current smoothness energy
    Eigen::VectorXd smoothness_sym;

    // Individual components 
    Eigen::VectorXd smoothness_L2;
    Eigen::VectorXd smoothness_L4;
    // Eigen::VectorXd smoothness_L2x2;

    // Individual components 
    Eigen::VectorXd curl_L2;
    Eigen::VectorXd curl_L4;
    // Eigen::VectorXd curl_L2x2;




};

enum class VariableType {
    primals, moments, deltas, gammas, Element_COUNT
};

// AppState holds the state of the application
class AppState {
public:

    ~AppState() {
        // delete config;
        delete os;
    }

    // File IO state 
    std::string directoryPath;  // this is where files are loaded from 
    std::string logFolderPath;  // this is where the output gets saved.  
                                // When load from a directory logFolderPath defaults to same directory, can change this.
    std::string meshName;
    // std::vector<std::string> bfraFiles;
    // std::vector<std::string> bmomFiles;
    std::optional<std::string> objFilePath;
    int currentFileID = -1;
    bool shouldReload = false; // this is a dynamic var which tells updaterendergeometry to reload data from the current directory
    bool updateRenderGeometryNextFrameIfPaused = false;

    // Init variables 
    Eigen::MatrixXd V; // Vertex positions
    Eigen::MatrixXi F; // Face indices

    /// @brief config optimization variables and params stored and initialized here.  ///
    std::unique_ptr<MyConfig> config;
    Eigen::VectorXi bound_face_idx;  
    // virtual Eigen::VectorXd dofs_to_vars(VariableType t); // TODO: implement this.
    // This would be a map that's used inside of OptZoo functions to make things more generic.  


    // Per solve metadata 
    std::string solveType = "SOLVER_TYPE_NOT_SET";
    std::string solveDescription = "SOLVER_DESCRIPTION_NOT_SET";

    std::string solveStatus = "STATUS_WAS_NOT_SET";
    std::string problemFileTag = "";



    // Optimization variables
    std::unique_ptr<Surface> cur_surf; // This initializes some more convenient data structures for building up local energies.
                      // In 3d need to use mesh data structures instead of surface.   
    bool keepSolving = true;
    int outerLoopIteration = 0;
    double cur_rel_residual = 0;
    double cur_abs_residual = 0;
    double cur_max_gradient_norm = 0;
    double cur_step_progress = 0;
    double cur_step_time = 0;

    // TODO: merge in the cube cover stuff.  
    Eigen::MatrixXd frames;
    Eigen::MatrixXd frames_orig;
    // Eigen::MatrixXd frames_viz;

    Eigen::MatrixXd moments; // TODO implement this!
    std::vector<Eigen::MatrixXd> frame_jacobians; // TODO implement this!
    Eigen::MatrixXd deltas;

    // Curl Operators 
    Eigen::MatrixXd C_primal;
    Eigen::MatrixXd C_sym_2;
    Eigen::MatrixXd C_sym_4; // TODO 

    // Per element selection indicies 
    // This is a much more generic way to create optzoo entries.  
    DOFMemoryLayout primals_layout;
    DOFMemoryLayout moments_layout;
    DOFMemoryLayout deltas_layout;

    // When you have already picked out the moments from the dofs 
    DOFMemoryLayout moments_L2_layout; // TODO implement support for this ! 
    DOFMemoryLayout moments_L4_layout;


    // Output State - gets visualized and serialized.  
    // could maybe also make this shared_ptr?
    OutputState* os; // TODO: make this a pointer to a base class, and then have derived classes for different output states.
                     // To makke this happen need to change the serialization code to support more generic file IO.  

    // GUI state 
    std::unordered_map<Field_View, FieldBounds> fieldBounds;
    // bool fieldViewActive [8] = {true, true, true, true, true, true, true, false};
    bool fieldViewActive [8] = {true, true, true, false,false,false, false, false};

    bool shouldLogData = true;
    Field_View prev_frame_element = Field_View::Element_COUNT;
    Field_View current_element;
    bool showVectorField = true;
    FieldBounds override_bounds; //= {0.0, 1e-4};
    bool override_bounds_active = false;
    bool show_frames = true;
    bool show_frames_as_lines = true;
    int max_saved_index = 0;

    // double L4_alpha = 0;
    float gui_vec_size = .01;
    Views::Sym_Moment_View cur_moment_view = Views::Sym_Moment_View::L2_plus_L4; 
    Views::Sym_Curl_View cur_curl_view = Views::Sym_Curl_View::L2_plus_L4; 


    bool LogToFile(const std::string suffix); // Log based on fieldViewActive state


    // simulation metadata 
    int currentIteration; 
    int maxIterations = 99999; // move to config
    // int innerLoopIteration; // move this to app state...
    double convergenceEpsilon = 1e-12;
    double convergenceThreshold; // not used right now...
    // double identityWeight;
    bool isConverged;

    // give these better names

    // Constructor
    AppState();



    // Methods for managing AppState
    void refreshFileLists();
    void selectFile(const std::string& filename);
    void refreshData();
    void updateSliderValue(int value);
    void serializeData();
    void deserializeData();
    void zeroPassiveVars();

    // Additional methods as needed
};

#endif // APPSTATE_H


