
///////////////////////
/////   init simulation
///////////////////////

#include "Mint2DHook.h"
#include "AppState.h"
#include "FileParser.h"
#include "Serialization.h"
#include <filesystem>
#include <sstream>
#include <iostream>
#include <chrono>
#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include "ImGuiWidgets.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"


#include <igl/readOBJ.h>
#include <igl/on_boundary.h>
#include "date.h"
#include "MyConfig.h"
#include "ADWrapper/ADFuncRunner.h"
#include "Surface.h"



void Mint2DHook::drawGUI() {


// original gui also 

/* 

      ImGui::InputDouble("Smoothness Weight", &w_smooth);
        ImGui::InputDouble("S Perp Weight", &w_s_perp);
        ImGui::InputDouble("Curl Weight", &w_curl);
        ImGui::InputDouble("Bound Weight", &w_bound);


/// Whatever maybe make this a dropdown eventually 
// From line 556 of imgui demo: https://skia.googlesource.com/external/github.com/ocornut/imgui/+/refs/tags/v1.73/imgui_demo.cpp
            const char* element_names[Field_View::Element_COUNT] = { "Vector Norms", "Delta Norms", "Vector Dirichlet", "Symmetric Dirichlet", "Vector Curl", "Symmetric Curl", "free" };
            const char* current_element_name = (current_element >= 0 && current_element < Field_View::Element_COUNT) ? element_names[current_element] : "Unknown";
            ImGui::PushItemWidth(300);
            ImGui::SliderInt("Shading Mode", (int *) &current_element, 0, Element_COUNT - 1, current_element_name);
            ImGui::PopItemWidth();
          

*/



    ImGuiWidgets::ShowMainGUI(*appState);
}

void Mint2DHook::updateRenderGeometry() {
    // Update visualization data based on current state in appState

    // appState->renderData.frames.resize(appState->simulationData.frames.rows(), 3);
    // appState->renderData.frames << appState->simulationData.frames, Eigen::MatrixXd::Zero(appState->simulationData.frames.rows(), 1);

    // appState->renderData.deltas = appState->simulationData.deltas;
    // appState->renderData.vec_curl = appState->simulationData.curls_primal;
    // appState->renderData.sym_curl = appState->simulationData.curls_sym;
    // appState->renderData.frame_smoothness = appState->simulationData.smoothness_primal;
    // appState->renderData.moment_smoothness = appState->simulationData.smoothness_sym;


    // Update the vertex positions and other data in Polyscope
    // polyscope::getSurfaceMesh("c")->updateVertexPositions(appState->V);
    // polyscope::getSurfaceMesh("c")->updateFaceIndices(appState->F);

    




////
//////  Here we calculate derived quantities. 
/// 
    appState->norms_vec = appState->frames.rowwise().norm();
    appState->norms_delta = appState->deltas.rowwise().norm();
    


    std::cout << "update render geometry" << std::endl;


    renderState = new AppState(*appState);

    renderState->frames.resize(renderState->frames.rows(), 3);
    renderState->frames << appState->frames, Eigen::MatrixXd::Zero(appState->frames.rows(), 1);


    // Log data if necessary
    if (appState->shouldLogData) {
        // Serialize and save the frame data
        appState->LogToFile();
        appState->currentFileID++;

        // Additional logging for any other fields in AppState as needed
    }


    // Update all of the calculated quantities of metadata.  


    // Advance simulation statistics in gui. 





    // // Handle the gui_free case for the Field_View enum
    // if (appState->currentElement == Field_View::gui_free) {
    //     // Implementation for gui_free case
    //     if (appState->customVisualizationEnabled) {
    //         // Custom visualization logic
    //     }
    // }

    // Request a redraw in Polyscope to update the visualization
    polyscope::requestRedraw();
}


    // Need to fill out viewer for each of: Field_View { vec_dirch, moment_dirch, sym_curl_residual, primal_curl_residual,
    void Mint2DHook::renderRenderGeometry()
    {
		// polyscope::getSurfaceMesh("c")->updateVertexPositions(renderP);
        
        // polyscope::getSurfaceMesh("c")->centerBoundingBox();
        // polyscope::getSurfaceMesh("c")->resetTransform();

        // // polyscope::getSurfaceMesh("c")->addFaceScalarQuantity("vec_norms", appState->frame_norms)->setEnabled(true);
        // polyscope::getSurfaceMesh("c")->addFaceScalarQuantity("vec_norms", appState->curls_primal)->setEnabled(true);


        
        // polyscope::requestRedraw();   

        // std::cout << "render renderGeometry" << std::endl;

        // polyscope::getSurfaceMesh("c")->updateVertexPositions(appState->V);
    // polyscope::getSurfaceMesh("c")->updateFaceIndices(appState->F);

    // Depending on the current element view, render different quantities
    switch (renderState->current_element) {
        case Field_View::vec_norms:
        //  std::cout << "render norms_vec" << std::endl;
            polyscope::getSurfaceMesh("c")->addFaceScalarQuantity("Norms Vector", renderState->norms_vec)->setEnabled(true);
            break;
        case Field_View::delta_norms:
        //  std::cout << "render norms_delta" << std::endl;
            polyscope::getSurfaceMesh("c")->addFaceScalarQuantity("Norms Delta", renderState->norms_delta)->setEnabled(true);
            break;
        case Field_View::vec_dirch:
        //  std::cout << "render smoothness_primal" << std::endl;
            polyscope::getSurfaceMesh("c")->addFaceScalarQuantity("Dirichlet Primal", renderState->smoothness_primal)->setEnabled(true);
            break;
        case Field_View::moment_dirch:
        //  std::cout << "render smoothness_sym" << std::endl;
            polyscope::getSurfaceMesh("c")->addFaceScalarQuantity("Dirichlet Moment", renderState->smoothness_sym)->setEnabled(true);
            break;
        case Field_View::primal_curl_residual:
        //  std::cout << "render primal_curl" << std::endl;
            polyscope::getSurfaceMesh("c")->addFaceScalarQuantity("Curl Primal Residual", renderState->curls_primal)->setEnabled(true);
            break;
        case Field_View::sym_curl_residual:
        //  std::cout << "render sym_curl" << std::endl;
            polyscope::getSurfaceMesh("c")->addFaceScalarQuantity("Curl Symmetric Residual", renderState->curls_sym)->setEnabled(true);
            break;
        case Field_View::gui_free:
            // Implement logic for gui_free if required
            break;
        default:
            std::cerr << "Unknown Field_View option selected in AppState." << std::endl;
            break;
    }

        // Update other visualization properties based on AppState
        // Example: Vector field visualization
        if (appState->showVectorField) {
            //  std::cout << "show vector field" << std::endl;

            //  std::cout << "appState->frames.rows() " << appState->frames.rows() << std::endl;
            //  std::cout << "renderState->frames.rows() " << renderState->frames.rows() << std::endl;
            //  polyscope::getSurfaceMesh("c");


            auto vectorField = polyscope::getSurfaceMesh("c")->addFaceVectorQuantity("Vector Field", renderState->frames);
            vectorField->setVectorColor(glm::vec3(0.7, 0.7, 0.7));
            vectorField->setEnabled(true);
        }

        polyscope::requestRedraw();
              




    }

void Mint2DHook::pause() {
    PhysicsHook::pause();
    appState->keepSolving = true;
}
    // Pause the simulati




void Mint2DHook::initSimulation() {
    // Load mesh using igl::readOBJ
        // igl::readOBJ("/home/josh/Documents/mint_redux/geometry-processing-starter-kit/tools/shared/" + cur_mesh_name + ".obj", V, F);

    if (appState->directoryPath.empty()) {
        std::cerr << "No directory path provided. TODO FIX THIS." << std::endl;
        appState->directoryPath = "../../results/BLAH"; // switch to ../shared/mint2d_testsequence
        // return;
    }
    appState->keepSolving = true;

    FileParser fileParser(appState->directoryPath);

    Eigen::MatrixXd V; // Temporary storage for vertices
    Eigen::MatrixXi F; // Temporary storage for faces

    // Check if .bfra and .bmom files exist
    // bool bfraExists = fileParser.parseLargestFile((appState->frames), FileType::BFRA);
    // bool bmomExists = fileParser.parseLargestFile((appState->deltas), FileType::BMOM);

    bool bfraExists = false;
    bool bmomExists = false;

    // Load mesh and config
    if (bfraExists && bmomExists) {
        // Deserialize configuration from a file
        // if (!Serialization::deserializeConfig(appState->config, appState->directoryPath + "/config.json")) {
        //     std::cerr << "Failed to load config from " << appState->directoryPath + "/config.json" << std::endl;
        //     // Handle error, possibly set default config
        // }
    } else {
        // Load default mesh and set default config
        std::string default_path = std::string(SOURCE_PATH) + "/../shared/" + appState->meshName + ".obj";
        
        // std::string default_path = "/home/josh/Documents/mint_redux/gpgpt/tools/shared/" + cur_mesh_name + ".obj";

        std::cout << default_path << std::endl;
        if (!igl::readOBJ(appState->objFilePath.value_or(default_path), V, F)) {
            std::cerr << "Failed to load mesh from " << appState->objFilePath.value_or(default_path) << std::endl;
            return;
        }
        appState->cur_surf = new Surface(V, F);

        // Set default configuration
        appState->config = new MyConfig(); // Assign default values to the config

        // // Serialize default config to a file
        // if (!Serialization::serializeConfig(appState->config, appState->directoryPath + "/config.json")) {
        //     std::cerr << "Failed to save default config to " << appState->directoryPath + "/config.json" << std::endl;
        // }
    }

    // fieldViewActive = 

    std::cout << "V.rows() " << V.rows() << "F.rows()" << F.rows() << std::endl;
 
    // Set mesh data to AppState
    appState->V = V;
    appState->F = F;

    this->resetAppState();


    // Initialize other parameters and logging folder
    // initializeOtherParameters();
    initializeLogFolder();

    // Register mesh with Polyscope
    polyscope::registerSurfaceMesh("c", appState->V, appState->F);
    polyscope::getSurfaceMesh("c")->setEdgeWidth(0.6);
    polyscope::view::resetCameraToHomeView();
}


/* 

    bool Mint2DHook::simulateOneStep()
    {
        // auto func = func_wrapper->get();
        if (cur_iter < max_iters)
        {
            cur_iter++;
            inner_loop_iter++;


            auto t1 = std::chrono::high_resolution_clock::now();


            // auto [f, g, H_proj] = func.eval_with_hessian_proj(x);
            auto [f, g, H_proj] = func.eval_with_derivatives(x);
            TINYAD_DEBUG_OUT("Energy in iteration " << cur_iter << ": " << f);
            // std::cout<<"the number of nonzeros "<<H_proj.nonZeros() << "number of non-zeros per dof " << H_proj.nonZeros() / (6*F.rows()) << " # rows " << H_proj.rows() << " faces " << F.rows() <<std::endl;

            // std::cout<<"the number of nonzeros "<<H_proj.nonZeros()<<std::endl;
            Eigen::VectorXd d;
            double dec;
            // d = TinyAD::newton_direction(g, H_proj, solver);
             // = TinyAD::newton_decrement(d, g);

            if (prev_energy < 0)
            {
              prev_energy = f + 100 * convergence_eps;
            }
            
            try
            {
              if (w_smooth_vector > 0 || useProjHessian)
              {
                auto [f_h, g_h, H_proj_h] = func.eval_with_hessian_proj(x);
                f = f_h;
                g = g_h;
                H_proj = H_proj_h;
                d = TinyAD::newton_direction(g, H_proj, solver, 0.);
                dec = TinyAD::newton_decrement(d, g);

                if ( dec / f < 1e-3)
                {
                  useProjHessian = false;
                  std::cout << "switch off projected hessian to fine-tune result" << std::endl;
                }


              }
              else
              {
                d = TinyAD::newton_direction(g, H_proj, solver, identity_weight);
                dec = TinyAD::newton_decrement(d, g);
                identity_weight = identity_weight / 2.;
              }
              
            }
            catch(const std::exception& e)
            {
              auto [f_h, g_h, H_proj_h] = func.eval_with_hessian_proj(x);
              f = f_h;
              g = g_h;
              H_proj = H_proj_h;
              d = TinyAD::newton_direction(g, H_proj, solver);
              dec = TinyAD::newton_decrement(d, g);
              if ( !useProjHessian )
                identity_weight = identity_weight * 10.;
            }
            
            // 
            std::cout << "current decrement: " << dec << std::endl;
            // if( dec < convergence_eps )
            // {
            //   buffer -= 1;
            //   identity_weight = identity_weight * 10.;
            // }

            if ( dec < convergence_eps || (inner_loop_iter > 300 && dec / f < 1e-5))
            {
              std::cout << "***** current decrement: " << dec << std::endl;
              buffer = 5;
              identity_weight = 1e-6;
              if (w_smooth_vector > 0)
              {
                w_smooth_vector = 0;
                
              }
              else {
                w_attenuate = w_attenuate / 10.;
                std::cout << "New attenuation value is set to: " << w_attenuate << std::endl;
                // inner_loop_iter = 0;
                if (w_attenuate < 1e-12)
                  cur_iter = max_iters;
              }

              useProjHessian = true;
                 

                // Eigen::MatrixXd tmp =TinyAD::to_passive(H_proj);
                //  igl::writeDMAT("converged_hessian.dmat",tmp,true);
            }

            // Eigen::MatrixXd tmp =TinyAD::to_passive(H_proj);
            // igl::writeDMAT("curr_hessian.dmat",tmp,true);

            // cur_iter = max_iters; // break
            x = TinyAD::line_search(x, d, f, g, func, 1., .8, 512, 1e-3);

            auto t2 = std::chrono::high_resolution_clock::now();
            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            std::cout << ms_int.count() << "ms\n";



            ///// Move this out 
            func.x_to_data(x, [&] (int f_idx, const Eigen::VectorXd& v) {
                frames.row(f_idx) = v.head<2>();
                // metadata.row(f_idx) = v.segment(2, 2);
                deltas.row(f_idx) = v.tail<4>();
                // if (bound_face_idx(f_idx) == 1)
                // {
                //   frames.row(f_idx) = frames_orig.row(f_idx);
                // }
                });

    


        }
        else if (cur_iter == max_iters) 
        {
            TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));
            cur_iter++;

             // FINAL LOGGING.  
        }
        else{
            this->pause();
        }


        
        return false;
    }



*/


bool Mint2DHook::simulateOneStep() {
    int max_iters = appState->maxIterations;
    int cur_iter = appState->currentIteration;
    double convergence_eps = appState->convergenceEpsilon;
    double cur_obj = opt->eval_func_at(opt->get_current_x());
    if (cur_iter < max_iters && cur_obj > convergence_eps && appState->keepSolving)
    {
        cur_iter++;
        std::cout << "*** GLOBAL STEP: " << cur_iter << "*** "  << std::endl;
        std::cout << "DOFS in opt" << opt->_cur_x.rows() << std::endl;
        std::cout << "nvars in opt" << opt->get_num_vars() << std::endl; 
        std::cout << "cur_obj: " <<  cur_obj  << " convergence_eps: " << convergence_eps << std::endl;

        opt->take_newton_step(opt->get_current_x());
        if (opt->_dec < convergence_eps)
        {
            std::cout << "**** Converged, early exit ****" << std::endl;
            std::cout << "Current Objective is " << opt->get_fval_at_x() << std::endl;

            appState->keepSolving = false;

            // In principle should load next config because an optimziation might have a schedule.
            //  cur_iter = max_iters; // break
            //  adhook.reset_params();
        }

        updateAppStateFromOptState();

            

     }
        else if (cur_iter == max_iters) 
        {
            TINYAD_DEBUG_OUT("Final energy: " << func.eval(opt->get_current_x()));
            cur_iter++;
        }
        else{
            std::cout << "**** Pause Simulation ****" << std::endl;
            this->pause();
        }

    appState->currentIteration = cur_iter;

    return false;


}


void Mint2DHook::resetAppState() {
    // Resetting simulation parameters to default or initial values
    appState->currentIteration = 0;
    appState->maxIterations = 5000; // Default maximum iterations
    appState->convergenceEpsilon = 1e-10;



    // Resetting mesh data
    appState->frames.setZero(appState->F.rows(), 2);
    appState->deltas.setZero(appState->F.rows(), 4);
    appState->moments.setZero(appState->F.rows(), 0);

    // Resetting derived quantities
    appState->norms_vec.setZero(appState->F.rows());
    appState->norms_delta.setZero(appState->F.rows());
    appState->curls_primal.setZero(appState->F.rows());
    appState->curls_sym.setZero(appState->F.rows());
    appState->smoothness_primal.setZero(appState->F.rows());
    appState->smoothness_sym.setZero(appState->F.rows());

    // Reinitialize the boundary conditions if needed
    initBoundaryConditions();

    // Optionally, re-register mesh with Polyscope if visualization needs a reset
    polyscope::removeAllStructures();
    polyscope::registerSurfaceMesh("c", appState->V, appState->F);
    polyscope::getSurfaceMesh("c")->setEdgeWidth(0.6);
}

    //     appState->renderFrames = Eigen::MatrixXd::Zero(F.rows(), 3);

    // appState->frames = Eigen::MatrixXd::Zero(F.rows(), 2);
    // appState->moments = Eigen::MatrixXd::Zero(F.rows(), 0);
    // appState->deltas = Eigen::MatrixXd::Zero(F.rows(), 4);
    // //   metadata = Eigen::MatrixXd::Zero(F.rows(), 2);
    // appState->norms_vec = Eigen::VectorXd::Zero(F.rows());
    // appState->norms_delta = Eigen::VectorXd::Zero(F.rows());
    // appState->curls_sym = Eigen::VectorXd::Zero(F.rows());
    // appState->curls_primal = Eigen::VectorXd::Zero(F.rows());
    // appState->smoothness_primal = Eigen::VectorXd::Zero(F.rows());
    // appState->smoothness_sym = Eigen::VectorXd::Zero(F.rows());



void Mint2DHook::initializeLogFolder() {
    // Retrieve current time
    auto now = std::chrono::system_clock::now();
    auto date = date::floor<date::days>(now);
    auto ymd = date::year_month_day{date};
    auto time = date::make_time(std::chrono::duration_cast<std::chrono::milliseconds>(now - date));

    // Construct the log folder path using the mesh name and the current date-time
    std::ostringstream folderStream;
    folderStream << "../../results/" << appState->meshName << "_"
                 << static_cast<unsigned>(ymd.month()) << "_"
                 << static_cast<unsigned>(ymd.day()) << "_"
                 << time.hours().count() << "_"
                 << time.minutes().count();

    // Store the constructed path in AppState
    appState->logFolderPath = folderStream.str();

    // Create the directory using std::filesystem
    std::filesystem::create_directories(appState->logFolderPath);

    // Inform the user about the log folder creation
    std::cout << "Log folder created at: " << appState->logFolderPath << std::endl;
}


void Mint2DHook::initializeOtherParameters() {
    // ... implementation for initializing other parameters ...
}



void Mint2DHook::initBoundaryConditions() {
    // Assuming boundary faces are identified in AppState
    Eigen::MatrixXi K;

    // Eigen::MatrixXi bound_face_idx = appState->bound_face_idx;

    Eigen::VectorXi boundaryFaces;
    igl::on_boundary(appState->F,boundaryFaces, K);

    appState->bound_face_idx = boundaryFaces;

    // Initialize boundary conditions
    for (int i = 0; i < boundaryFaces.size(); ++i) {
        if (boundaryFaces(i) == 1) { // If face is on the boundary
            Eigen::RowVector3d centroid = (appState->V.row(appState->F(i, 0)) +
                                           appState->V.row(appState->F(i, 1)) +
                                           appState->V.row(appState->F(i, 2))) / 3.0;

            if (centroid.norm() < 0.45) { // Custom condition for boundary faces
                boundaryFaces(i) = -1; // Mark for special handling or exclusion
            } else {
                // Set frame orientation based on the centroid
                Eigen::Vector2d frame = Eigen::Vector2d(centroid.y(), -centroid.x()).normalized();
                appState->frames.row(i) = frame;
            }
        }
    }

    appState->frames_orig = appState->frames;

    // 

    // 

}






void Mint2DHook::updateOptimizationParameters() {
    // Example implementation (needs to be adapted to specific needs)
    if (appState->currentIteration % 10 == 0) {
        // Adjust weights or parameters every 10 iterations
        appState->identityWeight *= 0.9; // Example: gradually decrease the identity weight
    }
    // Add other conditional parameter adjustments as needed
}


void Mint2DHook::checkAndUpdateConvergence(double decrement, double energy) {
    // Example convergence check
    if (decrement < appState->convergenceThreshold) {
        appState->isConverged = true;
        std::cout << "Optimization has converged." << std::endl;
    }
}

void Mint2DHook::finalizeIteration() {
    // Log current state if necessary
    // Check for additional stopping conditions
    // Example:
    if (appState->currentIteration >= appState->maxIterations) {
        std::cout << "Reached maximum number of iterations." << std::endl;
        this->pause(); // Pause the simulation
    }
}



/*
*
*
*
*

  bool boundaryConditionsLoaded = false;
    // Logic to check and load boundary conditions from a file
    // If successful:
    // boundaryConditionsLoaded = true;

    if (!boundaryConditionsLoaded) {
        // Set default boundary conditions
        Eigen::MatrixXd surfNormals; // Initialize with surface normals of your mesh
        Eigen::MatrixXd surfCenters; // Initialize with surface centers of your mesh

        // Call the radialBoundConstraints function
        Eigen::SparseMatrix<double> B_op_sparse;
        Eigen::VectorXd targ_moms;
        radialBoundConstraints(surfNormals, surfCenters, appState.nAugFrames, appState.logFolderPath + "/default_boundary.targ",
                               B_op_sparse, targ_moms);

        // Now B_op_sparse and targ_moms have default boundary conditions
        // Apply these conditions to your simulation as needed
    }



    void radialBoundConstraints(const Eigen::MatrixXd& surfNormals, const Eigen::MatrixXd& surfCenters, 
                            int nAugFrames, const std::string& logname,
                            Eigen::SparseMatrix<double>& B_op_sparse, Eigen::VectorXd& targ_moms);


void radialBoundConstraints(const Eigen::MatrixXd& surfNormals, const Eigen::MatrixXd& surfCenters, 
                            int nAugFrames, const std::string& logname,
                            Eigen::SparseMatrix<double>& B_op_sparse, Eigen::VectorXd& targ_moms) {
    // ... existing initialization logic ...

    // Processing surface normals and centers
    Eigen::MatrixXd selNormals = surfNormals;
    selNormals.col(2).setZero();
    Eigen::VectorXd normValues = selNormals.rowwise().norm();
    std::vector<int> curPinnedBoundIndices;

    for (int i = 0; i < normValues.size(); ++i) {
        if (normValues(i) > std::sqrt(0.001)) {
            curPinnedBoundIndices.push_back(i);
        }
    }

    Eigen::MatrixXd selCenters(curPinnedBoundIndices.size(), 3);
    for (size_t i = 0; i < curPinnedBoundIndices.size(); ++i) {
        selCenters.row(i) = surfCenters.row(curPinnedBoundIndices[i]);
    }

    int ntets = nAugFrames - static_cast<int>(selNormals.rows());
    Eigen::MatrixXd pinnedFrames(3, 3 * curPinnedBoundIndices.size());

    // Processing pinned frames
    for (size_t i = 0; i < curPinnedBoundIndices.size(); ++i) {
        Eigen::Vector3d tc = selCenters.row(i);
        double theta = std::atan2(tc.y(), tc.x());
        double r = tc.x() * tc.x() + tc.y() * tc.y();
        r = std::sqrt(r);

        // Rotation matrix
        Eigen::Matrix3d rot;
        rot << std::cos(theta), -std::sin(theta), 0,
               std::sin(theta), std::cos(theta), 0,
               0, 0, 1;

        // Apply transformation
        Eigen::Matrix3d f = rot; // Add your transformation logic here
        pinnedFrames.block<3, 3>(0, 3 * i) = f;
    }

    // TODO: Export to TARG file using `writeTARG` function and import using `readTARG` function

    // You need to implement writeTARG and readTARG functions to handle file IO
    // Example: writeTARG(logname, pinnedFrames);
    // Example: readTARG(logname, B_op_sparse, targ_moms);

    // Visualization logic if needed
    // visualizeFrameField(pinnedFrames, selCenters); // Implement as needed
}




*/