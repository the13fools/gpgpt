#include "AppState.h"
#include "FileParser.h"
#include "Serialization.h"
#include "FieldView.h"

#include <iostream>
#include <Eigen/Sparse>

#include <fstream>

// Constructor here
AppState::AppState()
{
    os = new OutputState();

    override_bounds.lower = 1e-16;
    override_bounds.upper = 1e-5;

    std::vector< Eigen::Triplet<double> > triplets;

    std::vector<double> L2_sym_weights = {1., 2., 2., 1., 2., 1.};
    std::vector<double> L4_sym_weights = { 1., 4., 4., 6., 12., 6., 4., 12., 12., 4., 1., 4., 6., 4., 1.};

    std::vector<double> L2_curl_weights = {1., 2., 1.};
    std::vector<double> L4_curl_weights =  {1., 4., 6., 4., 1.};
    std::vector<double> L6_curl_weights =  {1., 6., 15., 20., 15., 6., 1.};


    setSparseMetricFromWeights(L2_sym_tensor_weights, L2_sym_weights);
    setSparseMetricFromWeights(L4_sym_tensor_weights, L4_sym_weights);
    setSparseMetricFromWeights(L2_curl_tensor_weights, L2_curl_weights);
    setSparseMetricFromWeights(L4_curl_tensor_weights, L4_curl_weights);
    setSparseMetricFromWeights(L6_curl_tensor_weights, L6_curl_weights);
    // for (int i = 0; i < 6; i++)
    // {
    //     triplets.push_back(Eigen::Triplet<double>(i, i, L2_sym_weights.at(i)));
    // }

    // L2_sym_tensor_weights = Eigen::SparseMatrix<double>(6, 6);
    // L2_sym_tensor_weights.setFromTriplets(triplets.begin(), triplets.end());
    // L2_sym_tensor_weights.makeCompressed();

    // triplets.clear();

    // for (int i = 0; i < 15; i++)
    // {
    //     triplets.push_back(Eigen::Triplet<double>(i, i, L4_sym_weights.at(i)));
    // }

    // L4_sym_tensor_weights = Eigen::SparseMatrix<double>(15, 15);
    // L4_sym_tensor_weights.setFromTriplets(triplets.begin(), triplets.end());
    // L4_sym_tensor_weights.makeCompressed();

}

void AppState::setSparseMetricFromWeights(Eigen::SparseMatrix<double> &M, const std::vector<double> weights)
{
    std::vector< Eigen::Triplet<double> > triplets;

    for (int i = 0; i < weights.size(); i++)
    {
        triplets.push_back(Eigen::Triplet<double>(i, i, weights.at(i)));
    }

    M = Eigen::SparseMatrix<double>(weights.size(), weights.size());
    M.setFromTriplets(triplets.begin(), triplets.end());
    M.makeCompressed();
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
        
    int num_vecs = this->frames.size();

    std::cout << "num_vecs: " << num_vecs << std::endl;

// TODO add moments here 
    Serialization::serializeMatrix(this->frames, this->logFolderPath + "/frames"+ "_" + suffix + ".bfra", num_vecs);

    // TODO: Save frames prev and visualize diff
    // Serialization::serializeMatrix(this->deltas, this->logFolderPath + "/deltas" + "_" + suffix + ".bmom", num_vecs);
    // Serialization::serializeMatrix(this->deltas, this->logFolderPath + "/deltas" + "_" + suffix + ".bmom", num_vecs);



    Serialization::serializeConfig(*this->config, this->logFolderPath + "/config" + "_" + suffix + ".json");

    LogCurrentOptStats();

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

void AppState::LogCurrentOptStats()
{
    // std::ofstream out;
    // out.open("myfile.txt", std::ios::app);
    std::ofstream energyOut(this->logFolderPath + "/total_energy.txt", std::ios::app);
    energyOut << this->os->cur_global_objective_val << std::endl;
    energyOut.close();

    std::ofstream gradNormOut(this->logFolderPath + "/max_grad_norm.txt", std::ios::app);
    gradNormOut << std::scientific << this->cur_max_gradient_norm << std::endl;
    gradNormOut.close();

    // ofstream energyOut(this->logFolderPath + "/total_energy.txt", ios::app);
    // energyOut << this.cur_max_gradient_norm() << std::endl;
    // energyOut.close();
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
