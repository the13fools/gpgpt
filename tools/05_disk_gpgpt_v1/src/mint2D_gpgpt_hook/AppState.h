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
    Eigen::MatrixXd moments; // TODO implement this!
    Eigen::MatrixXd deltas;

    Eigen::VectorXd norms_vec;
    Eigen::VectorXd norms_delta;
    Eigen::VectorXd norms_moment;
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

    bool logToFile();

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

/*

       Serialization::serializeMatrix(appState->frames, appState->logFolderPath + "/frames.bfra");
        Serialization::serializeMatrix(appState->deltas, appState->logFolderPath + "/deltas.bmom");

    

          // Here, we'll also log relevant data to files based on the fieldViewActive flags

        // Log other Eigen::Vectors based on fieldViewActive flags
        for (int i = 0; i < Field_View::Element_COUNT; ++i) {
            if (appState->fieldViewActive[i]) {
                // Determine the file path based on the field view
                Field_View cfv = static_cast<Field_View>(i);
                std::string stub = fieldViewToFileStub(cvf);
                std::string filePath = appState->logFolderPath + "/" + stub + "_" + std::to_string(appState->currentFileID + 100000) + ".bfra"; // better to call these bdat

                // Serialize the corresponding Eigen::Vector
                switch (static_cast<Field_View>(i)) {
                    case Field_View::vec_norms:
                        Serialization::serializeVector(appState->norms_vec, filePath);
                        break;
                    case Field_View::delta_norms:
                        Serialization::serializeVector(appState->norms_delta, filePath);
                        break;
                    case Field_View::vec_dirch:
                        Serialization::serializeVector(appState->smoothness_primal, filePath);
                        break;
                    case Field_View::moment_dirch:
                        Serialization::serializeVector(appState->smoothness_sym, filePath);
                        break;
                    case Field_View::primal_curl_residual:
                        Serialization::serializeVector(appState->curls_primal, filePath);
                        break;
                    case Field_View::sym_curl_residual:
                        Serialization::serializeVector(appState->curls_sym, filePath);
                        break;
                    // ... handle other Field_View cases as needed
                    default:
                        break; // Unknown or unsupported field view
                }
            }
        }

*/
