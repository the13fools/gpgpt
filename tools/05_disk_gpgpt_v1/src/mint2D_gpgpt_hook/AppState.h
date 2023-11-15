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

// Enum for identifying field view quantities
enum Field_View {
    vec_norms, delta_norms, vec_dirch, moment_dirch, primal_curl_residual, sym_curl_residual, gui_free, Element_COUNT
};

// Struct to hold bounds for each field view quantity
struct FieldBounds {
    double upper = std::numeric_limits<double>::max();
    double lower = std::numeric_limits<double>::lowest();
};

// AppState holds the state of the application
class AppState {
public:
    std::string logFolderPath;
    std::string directoryPath;
    std::string meshName;
    std::vector<std::string> bfraFiles;
    std::vector<std::string> bmomFiles;
    std::optional<std::string> objFilePath;
    int currentFileID = 0;
    std::unordered_map<Field_View, FieldBounds> fieldBounds;
    bool fieldViewActive [8] = {false};
    bool shouldLogData = true;
    Field_View current_element;

    Eigen::VectorXi boundaryFaces;
    // Optimization variables
    Eigen::MatrixXd frames;
    Eigen::MatrixXd deltas;
    Eigen::VectorXd curls_primal;
    Eigen::VectorXd curls_sym;
    Eigen::VectorXd smoothness_primal;
    Eigen::VectorXd smoothness_sym;
    Eigen::MatrixXd V; // Vertex positions
    Eigen::MatrixXi F; // Face indices

    MyConfig* config;

// give these better names
    int currentIteration; 
    int maxIterations;
    int innerLoopIteration;
    double convergenceEpsilon;
    double convergenceThreshold;
    double identityWeight;
    bool isConverged;

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
