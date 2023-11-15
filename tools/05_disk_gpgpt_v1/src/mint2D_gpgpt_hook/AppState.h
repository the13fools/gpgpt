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

using Field_View = Views::Field_View;

// Enum for identifying field view quantities
// enum Field_View {
//     vec_norms, delta_norms, vec_dirch, moment_dirch, primal_curl_residual, sym_curl_residual, gui_free, Element_COUNT
// };
// Views::

// enum Field_View : unsigned int;

// Struct to hold bounds for each field view quantity
struct FieldBounds {
    float upper = std::numeric_limits<float>::max();
    float lower = std::numeric_limits<float>::lowest();
};

// AppState holds the state of the application
class AppState {
public:

    // File IO state 
    std::string directoryPath;  // this is where files are loaded from 
    std::string logFolderPath;  // this is where the output gets saved.  
                                // When load from a directory logFolderPath defaults to same directory, can change this.
    std::string meshName;
    std::vector<std::string> bfraFiles;
    std::vector<std::string> bmomFiles;
    std::optional<std::string> objFilePath;
    int currentFileID = 0;



    // Init variables 
    Eigen::MatrixXd V; // Vertex positions
    Eigen::MatrixXi F; // Face indices
    MyConfig* config;
    Eigen::VectorXi boundaryFaces;


    // Optimization variables
    Eigen::MatrixXd frames;
    Eigen::MatrixXd moments; // TODO implement this!
    Eigen::MatrixXd deltas;


    // derived variables 
    Eigen::VectorXd norms_vec;
    Eigen::VectorXd norms_delta;
    Eigen::VectorXd norms_moment; // TODO 
    Eigen::VectorXd curls_primal;
    Eigen::VectorXd curls_sym;
    Eigen::VectorXd smoothness_primal;
    Eigen::VectorXd smoothness_sym;


    // GUI state 
    std::unordered_map<Field_View, FieldBounds> fieldBounds;
    bool fieldViewActive [8] = {false};
    bool shouldLogData = true;
    Field_View current_element;
    bool showVectorField;

    bool LogToFile(); // Log based on fieldViewActive state


    // simulation metadata 
    int currentIteration; 
    int maxIterations;
    int innerLoopIteration;
    double convergenceEpsilon;
    double convergenceThreshold;
    double identityWeight;
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

    // Additional methods as needed
};

#endif // APPSTATE_H


