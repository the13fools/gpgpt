#include "AppState.h"
#include "FileParser.h"
#include "Serialization.h"
#include "FieldView.h"

#include <iostream>

// Constructor here
AppState::AppState()
{
    os = new OutputState();

    override_bounds.lower = 1e-16;
    override_bounds.upper = 1e-5;
}



// Implement the refreshFileLists to populate the lists of files
void AppState::refreshFileLists() {
    FileParser fileParser(directoryPath);



    // // Refresh bfraFiles and bmomFiles
    // fileParser.parseBfraFiles(bfraFiles);
    // fileParser.parseBmomFiles(bmomFiles);

    // Refresh objFilePath
    // fileParser.parseObjFile(objFilePath);

    // Additional logic to handle other types of files if needed
}


bool AppState::LogToFile(const std::string suffix)
{
        
    int num_vecs = this->os->frames.size();

    std::cout << "num_vecs: " << num_vecs << std::endl;

// TODO add moments here 
    Serialization::serializeMatrix(this->frames, this->logFolderPath + "/frames"+ "_" + suffix + ".bfra", num_vecs);
    Serialization::serializeMatrix(this->deltas, this->logFolderPath + "/deltas" + "_" + suffix + ".bmom", num_vecs);

    Serialization::serializeConfig(*this->config, this->logFolderPath + "/config" + "_" + suffix + ".json");


          // Here, we'll also log relevant data to files based on the fieldViewActive flags

        // Log other Eigen::Vectors based on fieldViewActive flags
        for (int i = 0; i < (int) Field_View::Element_COUNT; ++i) {
            if (this->fieldViewActive[i]) {
                // Determine the file path based on the field view
                Field_View cfv = static_cast<Field_View>(i);
                std::string stub = fieldViewToFileStub(cfv);
                std::string filePath = this->logFolderPath + "/" + stub + "_" + suffix + ".bdat"; // better to call these bdat

                // Serialize the corresponding Eigen::Vector
                switch (static_cast<Field_View>(i)) {
                    case Field_View::vec_norms:
                        Serialization::serializeVector(os->norms_vec, filePath);
                        break;
                    case Field_View::delta_norms:
                        Serialization::serializeVector(os->norms_delta, filePath);
                        break;
                    case Field_View::vec_dirch:
                        Serialization::serializeVector(os->smoothness_primal, filePath);
                        break;
                    case Field_View::moment_dirch:
                        // Serialization::serializeVector(os->smoothness_L2, filePath);
                        // TODO FIX THIS LOGGING
                        Serialization::serializeVector(os->smoothness_sym, filePath);
                        break;
                    case Field_View::primal_curl_residual:
                        Serialization::serializeVector(os->curls_primal, filePath);
                        break;
                    case Field_View::sym_curl_residual:
                        Serialization::serializeVector(os->curls_sym, filePath);
                        break;
                    // ... handle other Field_View cases as needed
                    default:
                        break; // Unknown or unsupported field view
                }
            }
        }


    return true;
}

void AppState::zeroPassiveVars()
{
    // zero out the passive variables
    // this->os->norms_vec.setZero();
    // this->os->norms_delta.setZero();
    // this->os->norms_moment.setZero();
    this->os->thetas.setZero();
    this->os->curls_primal.setZero();
    this->os->curls_sym.setZero();
    this->os->smoothness_primal.setZero();
    this->os->smoothness_sym.setZero();
    this->os->smoothness_L2.setZero();
    this->os->smoothness_L4.setZero();
    // this->os->smoothness_L2x2.setZero();
}

// // Method to log current variables (if logging is enabled)
// void AppState::logCurrentVariables() {
//     if (logVariables) {
//         // Implement serialization functions to log variables to a file
//         std::string logFilePath = directoryPath + "/log.json";

//         // // Serialize variables to log file
//         // if (Serialization::serializeAppState(*this, logFilePath)) {
//         //     std::cout << "Logged variables to: " << logFilePath << std::endl;
//         // } else {
//         //     std::cerr << "Failed to log variables." << std::endl;
//         // }
//     }
// }

// // Method to update the current frame and moment indices
// void AppState::updateIndices(int frameIndex, int momentIndex) {
//     currentFrameIndex = frameIndex;
//     currentMomentIndex = momentIndex;
//     // Additional logic to handle index updates
// }

// // Method to serialize the current app state to a file
// bool AppState::serializeData() {
//     std::string dataFilePath = directoryPath + "/appState.json";

//     // Serialize app state to a JSON file
//     if (Serialization::serializeAppState(*this, dataFilePath)) {
//         std::cout << "Serialized app state to: " << dataFilePath << std::endl;
//         return true;
//     } else {
//         std::cerr << "Failed to serialize app state." << std::endl;
//         return false;
//     }
// }

// // Method to deserialize app state from a file
// bool AppState::deserializeData() {
//     std::string dataFilePath = directoryPath + "/appState.json";

//     // Deserialize app state from a JSON file
//     if (Serialization::deserializeAppState(*this, dataFilePath)) {
//         std::cout << "Deserialized app state from: " << dataFilePath << std::endl;
//         return true;
//     } else {
//         std::cerr << "Failed to deserialize app state." << std::endl;
//         return false;
//     }
// }

// // Additional AppState-related functionality can be added here as needed.
