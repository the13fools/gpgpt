
///////////////////////
/////   init simulation
///////////////////////

#include "Mint3DHook.h"
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
#include "polyscope/volume_mesh.h"
#include "polyscope/point_cloud.h"


// #include <igl/readOBJ.h>
#include "CubeCover/readMesh.h"
#include <igl/writeOBJ.h>
#include <igl/on_boundary.h>
#include "date.h"
#include "MyConfig.h"
#include "ADWrapper/ADFuncRunner.h"
#include "Surface.h"
#include "UtilsMisc.h"

#include <thread>

void Mint3DHook::drawGUI() {


    // original gui also 

    /*

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

void Mint3DHook::updateRenderGeometry() {
    // Check if need to reload, and do this if needed.

    // Set the smoothness and curl to visualize according to the gui 
    switch (appState->cur_moment_view)
    {
        case Views::Sym_Moment_View::L2: 
            appState->os->smoothness_sym = appState->os->smoothness_L2;
            break;
        case Views::Sym_Moment_View::L4: 
            appState->os->smoothness_sym = appState->os->smoothness_L4;
            break;
        case Views::Sym_Moment_View::L2_plus_L4: 
            appState->os->smoothness_sym = appState->os->smoothness_L2 + appState->os->smoothness_L4;
            break;
    }

// This is kinda brittle... 
    int num_curls = appState->os->curls_Lks.size();
    switch (appState->cur_curl_view)
    {
        case Views::Sym_Curl_View::L2: 
            if(num_curls >= 1)
                appState->os->curls_sym = appState->os->curls_Lks[0];
            break;
        case Views::Sym_Curl_View::L4: 
            if(num_curls >= 2)
                appState->os->curls_sym = appState->os->curls_Lks[1];
            break;
        case Views::Sym_Curl_View::L6: 
            if(num_curls >= 3)
                appState->os->curls_sym = appState->os->curls_Lks[2];
            break;
        case Views::Sym_Curl_View::Total: 
            appState->os->curls_sym *= 0;
            for(int i = 0; i < num_curls; i++)
            {
                appState->os->curls_sym += appState->os->curls_Lks[i];
            }
            break;
    }


    if (appState->shouldReload) {
        std::cout << "do from file reload" << std::endl;
        pause();
        if (!this->isPaused())
        {
            std::cout << "failed to pause" << std::endl;
        }
        initSimulation();
    }
    appState->updateRenderGeometryNextFrameIfPaused = false;

    /// 
    appState->os->norms_vec = appState->frames.rowwise().norm();
    std::cout << "appState->os->norms_vec.rows(): " << appState->os->norms_vec.rows() <<  " maxcoeff: " << appState->os->norms_vec.maxCoeff() << std::endl;
    appState->os->norms_delta = appState->deltas.rowwise().norm();



    ////// The quantities not explicitly set above are set inside of to_passive calls in optzoo.
    // std::cout << opt->get_current_x() << std::endl;
    Eigen::VectorXd tmp = opt->get_current_x();
    std::cout << tmp.rows() << " " << tmp.cols() << std::endl;

  
    

    if ( tmp.rows() == 0 )
    {
        std::cout << "ERROR INITIALIZING OPTIMIZATION PROBLEM - OPT STATE CURRENTLY EMPTY" << std::endl;
        // tmp = Eigen::VectorXd::Zero(opt->get_num_vars());
    }
    else
    {
        try {

            //// Compute different parts of the objective function by changing hyper parameters to
            //// switch things on and off.  

            // A little fragile, make sure these are correct if you change stuff. 

            double prev_curl = appState->config->w_curl;
            double prev_smoothness = appState->config->w_smooth;
            double prev_attenuate = appState->config->w_attenuate;
            double prev_bound = appState->config->w_bound;

            appState->os->cur_global_objective_val = opt->eval_func_at(tmp);

            appState->config->w_curl = 0.;
            appState->config->w_smooth = 1.;
            appState->config->w_attenuate = 1.;
            appState->config->w_bound = 0;

            appState->obj_smoothness_part = opt->eval_func_at(tmp);

            appState->config->w_curl = 1.;
            appState->config->w_smooth = 0.;

            appState->obj_curl_part = opt->eval_func_at(tmp);

            appState->config->w_curl = prev_curl;
            appState->config->w_smooth = prev_smoothness;
            appState->config->w_attenuate = prev_attenuate;
            appState->config->w_bound = prev_bound;
            // appState->energy_trace.push_back(cur_obj);

            // opt->eval_func_with_derivatives(tmp);

            // appState->cur_abs_residual = opt->_dec;
            // appState->cur_max_gradient_norm = opt->_max_gradient_norm;

            // std::cout << "cur_obj " << appState->os->cur_global_objective_val << std::endl;
        }
        catch (std::runtime_error& e)
        {
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << "optimization crashing when evaluating the objective " << e.what() << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

        }

    }
    std::cout << "update render geometry" << std::endl;

    outputData = new OutputState(*appState->os);

    // clear previous passive_vars



    // renderState = new AppState(*appState);

    outputData->frames.clear();
    outputData->boundary_frames.clear();

    int frame_rank = appState->primals_layout.rank;
    int vec_dofs = appState->primals_layout.size / frame_rank;
    int frows = appState->frames.rows();
    int bfrows = appState->boundary_frames.rows();

    if (vec_dofs != 3 )
        std::cout << "******** wrong number of dofs for mint3d ********" << std::endl;

    // outputData->frames = appState->frames;

    for(int vid = 0; vid < frame_rank; vid++)
    {
        Eigen::MatrixXd vec_cur = Eigen::MatrixXd::Zero(appState->frames.rows(), 3);


        vec_cur << appState->frames.block(0, vid * vec_dofs, frows, vec_dofs); // ,  Eigen::MatrixXd::Zero(frows, 1);
        outputData->frames.push_back(vec_cur);
    }

    for(int vid = 0; vid < frame_rank; vid++)
    {
        Eigen::MatrixXd vec_cur = Eigen::MatrixXd::Zero(appState->boundary_frames.rows(), 3);
        vec_cur << appState->boundary_frames.block(0, vid * vec_dofs, bfrows, vec_dofs); // ,  Eigen::MatrixXd::Zero(frows, 1);
        outputData->boundary_frames.push_back(vec_cur);
    }

    // outputData->frames.resize(appState->frames.rows(), 3);
    // outputData->frames << appState->frames, Eigen::MatrixXd::Zero(appState->frames.rows(), 1);

    // if (appState->shouldReload)
    // {
    //     appState->currentFileID--;
    //     appState->shouldReload = false;
    // }

    // // Log data if necessary
    // if (appState->shouldLogData) {
    //     // Serialize and save the frame data
    //     // appState->currentFileID++;
    //     appState->max_saved_index = appState->max_saved_index + 1;
    //     std::string suffix = std::to_string(appState->max_saved_index + 100000);
    //     appState->LogToFile(suffix);

    //     if (appState->max_saved_index == appState->currentFileID + 1)
    //     {
    //         appState->currentFileID = appState->max_saved_index;
    //     }


    //     // Additional logging for any other fields in AppState as needed
    // }

     if (appState->shouldReload)
    {
        appState->currentFileID--;
        appState->shouldReload = false;
    }

 
     // Log data if necessary
     if (appState->shouldLogData) {
         // Serialize and save the frame data
        appState->currentFileID++;
        std::string suffix = std::to_string(appState->currentFileID + 100000);
         appState->LogToFile(suffix);
 
        if (appState->max_saved_index < appState->currentFileID)
        {
            appState->max_saved_index = appState->currentFileID;
        }

        appState->energy_trace.push_back(appState->os->cur_global_objective_val);
        appState->energy_smoothness_part_trace.push_back(appState->obj_smoothness_part);
        appState->energy_curl_part_trace.push_back(appState->obj_curl_part);
        appState->cur_max_gradient_norm_trace.push_back(appState->cur_max_gradient_norm);
        appState->solve_rel_residual_trace.push_back(appState->solve_rel_residual);
        appState->identity_weight_trace.push_back(appState->identity_weight);

     }


    appState->LogToFile("curr");


    // appState->zeroPassiveVars();

    // Request a redraw in Polyscope to update the visualization
    polyscope::requestRedraw();
}


// Need to fill out viewer for each of: Field_View { vec_dirch, moment_dirch, sym_curl_residual, primal_curl_residual,
void Mint3DHook::renderRenderGeometry()
{

    if ((appState->shouldReload || appState->updateRenderGeometryNextFrameIfPaused) && this->isPaused())
    {
        updateRenderGeometry();
        appState->updateRenderGeometryNextFrameIfPaused = false;
    }

    // Depending on the current element view, render different quantities
    const std::string cur_field = fieldViewToFileStub(appState->current_element);


    Eigen::VectorXd cur_scalar_quantity;
    switch (appState->current_element) {
        case Field_View::vec_norms:
            cur_scalar_quantity = outputData->norms_vec;
            break;
        case Field_View::delta_norms:
            cur_scalar_quantity = outputData->norms_delta;
            break;
        case Field_View::vec_dirch:
            cur_scalar_quantity = outputData->smoothness_primal;
            break;
        case Field_View::moment_dirch:
            // cur_scalar_quantity = outputData->smoothness_sym;
            // cur_scalar_quantity = outputData->smoothness_L2;
            cur_scalar_quantity = outputData->smoothness_sym;
            break;
        case Field_View::primal_curl_residual:
            cur_scalar_quantity = outputData->curls_primal;
            break;
        case Field_View::sym_curl_residual:
            cur_scalar_quantity = outputData->curls_sym;
            break;
        case Field_View::gui_free:
            // Implement logic for gui_free if required
            break;
        default:
            std::cerr << "Unknown Field_View option selected in AppState." << std::endl;
            break;
    }


    
    if (appState->current_element == Field_View::gui_free)
    {
        // noop
    }
    else
    {
        auto cur_scalar_field = polyscope::getVolumeMesh("c")->addCellScalarQuantity(cur_field, cur_scalar_quantity);
        if (appState->prev_frame_element != appState->current_element)
        {
            cur_scalar_field->setEnabled(true);
            double maxbound = cur_scalar_quantity.maxCoeff();
            double minbound = cur_scalar_quantity.minCoeff();
            double diff = maxbound - minbound;
            FieldBounds cur_bounds = appState->fieldBounds[appState->current_element];
            double abs_max = cur_bounds.upper * diff + minbound;
            double abs_min = cur_bounds.lower * diff + minbound;

            if (appState->override_bounds_active)
            {
                abs_max = appState->override_bounds.upper;
                abs_min = appState->override_bounds.lower;
            }

            cur_scalar_field->setMapRange(std::make_pair(abs_min, abs_max));


        }
    }


    if (appState->current_element == Field_View::gui_free)
    {
        // noop
    }
    else if (appState->showVectorField) {
        
        int num_vecs = outputData->frames.size();
        for(int v = 0; v < num_vecs; v++)
        {

            Eigen::MatrixXd cur_vec = outputData->frames[v];

            double color_shift = (v+1.) * 1.0 / num_vecs;

            auto vectorField = polyscope::getPointCloud("c_vecs")->addVectorQuantity("Vector Field " + std::to_string(v), cur_vec);
            auto vectorFieldNeg = polyscope::getPointCloud("c_vecs")->addVectorQuantity("Vector Field (negative) " + std::to_string(v), (-1.) * cur_vec);

            if(appState->show_frames && appState->show_frames_as_lines)
            {
                vectorField->setEnabled(true);
                vectorFieldNeg->setEnabled(true);
                vectorField->setVectorLengthScale(appState->gui_vec_size);
                vectorFieldNeg->setVectorLengthScale(appState->gui_vec_size);
                
            }
            else if (appState->show_frames)
            {
                vectorField->setEnabled(true);
                vectorFieldNeg->setEnabled(false);
            }
            else 
            {
                vectorField->setEnabled(false);
                vectorFieldNeg->setEnabled(false);
            }


        }

        for(int v = 0; v < num_vecs; v++)
        {

            Eigen::MatrixXd cur_vec = outputData->boundary_frames[v];

            double color_shift = (v+1.) * 1.0 / num_vecs;

            auto vectorField = polyscope::getPointCloud("surface_vecs")->addVectorQuantity("Vector Field " + std::to_string(v), cur_vec);
            auto vectorFieldNeg = polyscope::getPointCloud("surface_vecs")->addVectorQuantity("Vector Field (negative) " + std::to_string(v), (-1.) * cur_vec);

            if(appState->show_frames && appState->show_frames_as_lines)
            {
                vectorField->setEnabled(true);
                vectorFieldNeg->setEnabled(true);
                vectorField->setVectorLengthScale(appState->gui_vec_size);
                vectorFieldNeg->setVectorLengthScale(appState->gui_vec_size);
                
            }
            else if (appState->show_frames)
            {
                vectorField->setEnabled(true);
                vectorFieldNeg->setEnabled(false);
            }
            else 
            {
                vectorField->setEnabled(false);
                vectorFieldNeg->setEnabled(false);
            }


        }
    }

    

    polyscope::requestRedraw();





}

void Mint3DHook::pause() {
    PhysicsHook::pause();
    appState->keepSolving = true;
    appState->solveStatus = "paused";
    std::cout << "paused" << std::endl;
}
// Pause the simulation




void Mint3DHook::initSimulation() {
    // Load mesh using igl::readOBJ
        // igl::readOBJ("/home/josh/Documents/mint_redux/geometry-processing-starter-kit/tools/shared/" + cur_mesh_name + ".obj", V, F);

    appState->keepSolving = true;
    appState->solveStatus = "init simulation";


    bool create_new_dir = false;
    if (appState->directoryPath.empty()) {
        std::cerr << "No directory path provided. Creating new directory." << std::endl;
        create_new_dir = true;
        appState->directoryPath = "../../results/BLAH"; // switch to ../shared/mint2d_testsequence
        // return;
    }



    Eigen::MatrixXd V; // Temporary storage for vertices
    Eigen::MatrixXi T; // Temporary storage for tetrahedra
    Eigen::MatrixXi F; // Temporary storage for faces

    // Load mesh and config
    if (create_new_dir) {
        // Load default mesh and set default config
        std::string default_path = std::string(SOURCE_PATH) + "/../shared/" + appState->meshName + ".mesh";

        // std::string default_path = "/home/josh/Documents/mint_redux/gpgpt/tools/shared/" + cur_mesh_name + ".obj";

        std::cout << default_path << std::endl;
        
        if (!CubeCover::readMESH(appState->meshFilePath.value_or(default_path), V, T, F)) {
            std::cerr << "Failed to load mesh from " << appState->meshFilePath.value_or(default_path) << std::endl;

            // tricky part for windows
            default_path = std::string(SOURCE_PATH) + "/tools/shared/" + appState->meshName + ".mesh";
            if ( !CubeCover::readMESH(appState->meshFilePath.value_or(default_path), V, T, F) ) {
                std::cerr << "Failed to load mesh from " << appState->meshFilePath.value_or(default_path) << std::endl;
                return;
            }

        }
        appState->cur_surf = std::make_unique<Surface>(V, F);
        appState->cur_tet_mesh = std::make_unique<CubeCover::TetMeshConnectivity>(T);

        // Set default configuration
        appState->config = std::make_unique<MyConfig>(); // Assign default values to the config
        initConfigValues();
        // // Serialize default config to a file
        // if (!Serialization::serializeConfig(appState->config, appState->directoryPath + "/config.json")) {
        //     std::cerr << "Failed to save default config to " << appState->directoryPath + "/config.json" << std::endl;
        // }

    }
    else {

        // fileParser = new FileParser(appState->directoryPath);
        fileParser = std::make_unique<FileParser>(appState->directoryPath);
        appState->max_saved_index = fileParser->maxID;

        // TODO add parsing for this.
        // appState->config = new MyConfig(); // ADD FILE PARSING NOW! 

        // 

        /// TODO: FIX THIS
        std::string default_path = std::string(SOURCE_PATH) + "/../shared/" + appState->meshName + ".mesh";
        fileParser->meshFilePath = default_path; // 
        appState->meshFilePath = fileParser->meshFilePath;

        // TODO: Fix this to load from file.
        std::cout << "fileParser->meshFilePath " << fileParser->meshFilePath << std::endl;
        if ( !CubeCover::readMESH(fileParser->meshFilePath, V, T, F) ) {
            std::cerr << "Failed to load mesh from " << fileParser->meshFilePath << std::endl;
            return;
        }
        // appState->cur_surf = std::make_unique<Surface>(V, F);
        appState->cur_surf = std::make_unique<Surface>(V, F);
        appState->cur_tet_mesh = std::make_unique<CubeCover::TetMeshConnectivity>(T);

    }

    std::cout << "V.rows() " << V.rows() << " T.rows() " << T.rows() << " F.rows() " << F.rows() << std::endl;


    // Set mesh data to AppState
    appState->V = V;
    appState->T = T;
    appState->F = F;

    appState->nelem = T.rows();
    if ( appState->useBoundaryFrames )
    {
        appState->nelem += appState->cur_tet_mesh->nBoundaryElements();
    }

    if (appState->shouldReload)
    {

        if (appState->currentFileID == -1)
        {
            this->resetAppState();
            initializeLogFolder();
            appState->currentFileID = appState->max_saved_index;
            appState->config = std::make_unique<MyConfig>(); // Assign default values to the config
            initConfigValues();
        }

        

        loadPrimaryData();
        loadGuiState();

    }

    if (create_new_dir)
    {
        this->resetAppState();
        initializeLogFolder();
        appState->directoryPath = appState->logFolderPath;
        appState->currentFileID = 0;

    }

    // // save OBJ file 
    // if (!igl::writeOBJ(appState->directoryPath + "/" + appState->meshName + ".obj", V, F)) {
    //     std::cerr << "Failed to save mesh to " << appState->logFolderPath << std::endl;
    //     // return;
    // }

    // calculate tet centroids 
    // appState->tet_centroids = CubeCover::computeTetCentroids(appState->V, appState->T);

    // Register mesh with Polyscope
    // polyscope::registerTetMesh("c", appState->V, appState->T);






    // polyscope::view::resetCameraToHomeView();
}

bool Mint3DHook::loadPrimaryData() {
    bool success = true;

    std::vector<FileTypeInfo> fileInfo = {
        {"frames_", ".bfra", appState->frames},
        {"deltas_", ".bmom", appState->deltas},
        {"moments_", ".bmom", appState->moments}
    };

    for (const auto& info : fileInfo) {
        std::string file = fileParser->getFileWithID(info.prefix, info.extension, appState->currentFileID);
        // std::cout << file << std::endl;
        if (!file.empty()) {
            if (!Serialization::deserializeMatrix(info.targetMatrix, file)) {
                std::cerr << "Failed to load data for " << info.prefix << " from file: " << file << std::endl;
                success = false;
            }
        }
        else {
            std::cerr << "File not found for " << info.prefix << " with ID: " << appState->currentFileID << std::endl;
            success = false;
        }
    }

    std::string config_path = fileParser->getFileWithID("config_", ".json", appState->currentFileID);
    std::cout << "config_path " << config_path << std::endl;
    if (!config_path.empty()) {
        if (!Serialization::deserializeConfig(*appState->config, config_path)) {
            std::cerr << "Failed to load data for " << "config_" << " from file: " << config_path << std::endl;
            success = false;
        }
    }
    else {
        std::cerr << "File not found for " << "config_" << " with ID: " << appState->currentFileID << std::endl;
        success = false;
    }





    return success;
}



bool Mint3DHook::loadGuiState() {
    bool success = true;
    for (int i = 0; i < (int) Views::Field_View::gui_free; ++i) {
        Field_View view = static_cast<Field_View>(i);
        std::string fieldStub = Views::fieldViewToFileStub(view) + "_";
        std::string fieldFile = fileParser->getFileWithID(fieldStub, ".bdat", appState->currentFileID);

        if (!fieldFile.empty()) {
            Eigen::VectorXd data;
            if (!Serialization::deserializeVector(data, fieldFile)) {
                std::cerr << "Failed to load data for " << fieldStub << " from file: " << fieldFile << std::endl;
                success = false;
            }
            else {
                // Update the correct field in OutputState based on 'view'
                // Example: appState.os->norms_vec = data;
            }
        }
        else {
            std::cerr << "File not found for " << fieldStub << " with ID: " << appState->currentFileID << std::endl;
            success = false;
        }
    }
    return success;
}




bool Mint3DHook::simulateOneStep() {
    appState->solveStatus = "simulate one step";
    int max_iters = appState->maxIterations;
    int cur_iter = appState->currentIteration;
    double convergence_eps = appState->convergenceEpsilon;
    double cur_obj = opt->eval_func_at(opt->get_current_x());
    if (cur_iter < max_iters && cur_obj > convergence_eps && appState->keepSolving)
    {
        cur_iter++;
        std::cout << std::endl << "*** GLOBAL STEP: " << cur_iter << "*** " << std::endl;
        std::cout << "DOFS in opt" << opt->_cur_x.rows() << std::endl;
        std::cout << "nvars in opt" << opt->get_num_vars() << std::endl;

        // std::cout << "current x: " << opt->get_current_x().transpose() << std::endl;
        opt->take_newton_step(opt->get_current_x());
        double rel_res_correction = 1. / (cur_obj);
        appState->cur_rel_residual = std::max(opt->_dec * rel_res_correction, 1e-13);
        appState->cur_abs_residual = opt->_dec;
        appState->cur_max_gradient_norm = opt->_max_gradient_norm;
        appState->cur_step_progress = opt->_prev_step_progress;
        appState->cur_step_time = opt->_prev_step_time;
        appState->solve_residual = opt->_solve_residual;
        appState->rhs_norm = opt->_rhs_norm;
        appState->solve_rel_residual = opt->_solve_residual / opt->_rhs_norm;
        appState->identity_weight = opt->identity_weight;
        // std::cout << "cur_obj: " <<  cur_obj  << " conv_eps: " << convergence_eps << "rel_res_correction: " << rel_res_correction << std::endl;






        // Three conditions for convergence:
        // 1. Relative residual is small
        // 2. Absolute residual is small
        // 3. Gradient norm is small or step progress is vanishing/negative (line search failing to make progress)
        // if (appState->cur_rel_residual  < convergence_eps && appState->cur_abs_residual < 1e-4 && (appState->cur_max_gradient_norm < 1e-8 || opt->identity_vanished || opt->_prev_step_progress < 1e-10)) 
        if (appState->cur_max_gradient_norm < 1e-8 || opt->identity_weight > 1e12)
        {
            std::cout << "**** Converged current step ****" << std::endl;
            std::cout << "Current Objective is " << opt->get_fval_at_x() << std::endl;

            appState->keepSolving = false;

            // In principle should load next config because an optimziation might have a schedule.
            //  cur_iter = max_iters; // break
            //  adhook.reset_params();
            opt->reset_params();
        }

        updateAppStateFromOptState();
        appState->solveStatus = "finished one step";



    }
    else if (cur_iter == max_iters)
    {
        // TINYAD_DEBUG_OUT("Final energy: " << func.eval(opt->get_current_x()));
        cur_iter++;
    }
    else {
        std::cout << "**** Pause Simulation ****" << std::endl;
        this->pause();
    }

    appState->currentIteration = cur_iter;

    return false;


}

/// @brief This should only be called after the appstate has been initialized with a mesh 
void Mint3DHook::resetAppState() {
    // Resetting simulation parameters to default or initial values
    appState->currentIteration = 0;
    appState->max_saved_index = 0;

    // appState->currentFileID = 0;
    appState->maxIterations = 9999; // Default maximum iterations
    appState->convergenceEpsilon = 1e-12;// 1e-9;
    appState->outerLoopIteration = 0;

    appState->override_bounds.lower = 0;
    appState->override_bounds.upper = 1e-5;
    // appState->shouldReload = false;



    // Resetting mesh data
    appState->frames.setZero(appState->T.rows(), 3); // TODO make this generic 
    appState->boundary_frames.setZero(appState->cur_tet_mesh->nBoundaryElements(), 3); // TODO make this generic 

    appState->deltas.setZero(appState->T.rows(), 4);
    appState->moments.setZero(appState->T.rows(), 0);

    // Resetting derived quantities
    appState->os->norms_vec.setZero(appState->T.rows());
    appState->os->norms_delta.setZero(appState->T.rows());
    appState->os->curls_primal.setZero(appState->T.rows());
    appState->os->curls_sym.setZero(appState->T.rows());
    appState->os->smoothness_primal.setZero(appState->T.rows());
    appState->os->smoothness_L2.setZero(appState->T.rows());
    appState->os->smoothness_L4.setZero(appState->T.rows());
    appState->os->curl_L2.setZero(appState->T.rows());
    appState->os->curl_L4.setZero(appState->T.rows());
    // appState->os->smoothness_L2x2.setZero(appState->F.rows());

    // std::cout << "appState->os->norms_vec.rows(): " << appState->os->norms_vec.rows() << std::endl;

    appState->prev_frame_element = Field_View::Element_COUNT;

    // Reinitialize the boundary conditions if needed
    initBoundaryConditions();

    // Initialize symmetric curl operators 
    initCurlOperators();

    // Optionally, re-register mesh with Polyscope if visualization needs a reset
    polyscope::removeAllStructures();
    polyscope::registerTetMesh("c", appState->V, appState->T);
    polyscope::getVolumeMesh("c")->setEdgeWidth(0.6)->setTransparency(0.5);
    // polyscope::getVolumeMesh("c")->setEdgeWidth(0.6);
    polyscope::view::resetCameraToHomeView();


    appState->tet_centroids = Eigen::MatrixXd::Zero(appState->T.rows(), 3);
    for(int i = 0; i < appState->T.rows(); i++)
    {
        appState->tet_centroids.row(i) = (appState->V.row(appState->T(i, 0)) +
                                            appState->V.row(appState->T(i, 1)) +
                                            appState->V.row(appState->T(i, 2)) +
                                            appState->V.row(appState->T(i, 3))) / 4.0;
    }

    polyscope::registerPointCloud("c_vecs", appState->tet_centroids)->setPointRadius(0.0);


    int num_vecs = 3; // appState->frames.size();
    for(int v = 0; v < num_vecs; v++)
    {
        // Eigen::MatrixXd cur_vec = outputData->frames[v];
        Eigen::MatrixXd cur_vec = Eigen::MatrixXd::Zero(appState->frames.rows(), 3);

        double color_shift = (v+1.) * 1.0 / num_vecs;

        auto vectorField = polyscope::getPointCloud("c_vecs")->addVectorQuantity("Vector Field " + std::to_string(v), cur_vec);
            
        auto vectorFieldNeg = polyscope::getPointCloud("c_vecs")->addVectorQuantity("Vector Field (negative) " + std::to_string(v), (-1.) * cur_vec);
        vectorField->setVectorColor(glm::vec3(color_shift, 0.1, 0.3));
        vectorFieldNeg->setVectorColor(glm::vec3(color_shift, 0.9, 0.7));
        vectorField->setVectorRadius(0.01);
        vectorFieldNeg->setVectorRadius(0.01);
    }


    appState->bound_centroids = Eigen::MatrixXd::Zero(appState->cur_tet_mesh->nBoundaryElements(), 3);
    appState->bound_normals = appState->bound_centroids;
    appState->bound_b1 = appState->bound_centroids;
    appState->bound_b2 = appState->bound_centroids;
    for(int i = 0; i < appState->bound_centroids.rows(); i++)
    {
        int boundaryFace = appState->cur_tet_mesh->boundaryFace(i);
        Eigen::Vector3d a = appState->V.row(appState->cur_tet_mesh->faceVertex(boundaryFace, 0));
        Eigen::Vector3d b = appState->V.row(appState->cur_tet_mesh->faceVertex(boundaryFace, 1));
        Eigen::Vector3d c = appState->V.row(appState->cur_tet_mesh->faceVertex(boundaryFace, 2));

        Eigen::Vector3d b1 = (b - a).normalized();
        Eigen::Vector3d b2 = (c - a).normalized();
        Eigen::Vector3d n = b1.cross(b2).normalized();


        appState->bound_centroids.row(i) = ( a + b + c )  / 3.0;
        appState->bound_normals.row(i) = n;
        appState->bound_b1.row(i) = b1;
        appState->bound_b2.row(i) = n.cross(b1);
    }

    polyscope::registerPointCloud("surface_vecs", appState->bound_centroids)->setPointRadius(0.0);


    // int num_vecs = 3; // appState->frames.size();
    for(int v = 0; v < num_vecs; v++)
    {
        // Eigen::MatrixXd cur_vec = outputData->frames[v];
        Eigen::MatrixXd cur_vec = Eigen::MatrixXd::Zero(appState->boundary_frames.rows(), 3);

        double color_shift = (v+1.) * 1.0 / num_vecs;

        auto vectorField = polyscope::getPointCloud("surface_vecs")->addVectorQuantity("Vector Field " + std::to_string(v), cur_vec);
            
        auto vectorFieldNeg = polyscope::getPointCloud("surface_vecs")->addVectorQuantity("Vector Field (negative) " + std::to_string(v), (-1.) * cur_vec);
        vectorField->setVectorColor(glm::vec3(color_shift, 0.1, 0.3));
        vectorFieldNeg->setVectorColor(glm::vec3(color_shift, 0.9, 0.7));
        vectorField->setVectorRadius(0.01);
        vectorFieldNeg->setVectorRadius(0.01);
    }


    
}


void Mint3DHook::initCurlOperators()
{
    // this can probably be deleted but leaving in for now. 
    int ntets = appState->cur_tet_mesh->nTets();
    appState->R_facet_to_template.resize(ntets);
    appState->tet_facet_basis.resize(ntets);

    for (int tetidx = 0; tetidx < ntets; tetidx++)
    {
        std::vector<Eigen::MatrixXd> rot_facets_to_template;
        std::vector<Eigen::MatrixXd> cur_tet_bases;
        rot_facets_to_template.resize(4);
        cur_tet_bases.resize(4);
        for (int idx = 0; idx < 4; idx++)
        {
            // std::cout << "ntets " << appState->cur_tet_mesh->nTets() << "nfaces" << appState->cur_tet_mesh->nFaces() << std::endl;
            int face_idx = appState->cur_tet_mesh->tetFace(tetidx, idx);

            Eigen::Matrix3d rot_facet_to_template;
            Eigen::MatrixXd cur_facet_basis = Eigen::MatrixXd::Zero(2, 3);

            if (face_idx == -1)
            {
                rot_facet_to_template.setIdentity();
                rot_facets_to_template.at(idx) = rot_facet_to_template;
                continue;
            }

            Eigen::VectorXd a = appState->V.row(appState->cur_tet_mesh->faceVertex(face_idx, 0));
            Eigen::VectorXd b = appState->V.row(appState->cur_tet_mesh->faceVertex(face_idx, 1));
            Eigen::VectorXd c = appState->V.row(appState->cur_tet_mesh->faceVertex(face_idx, 2));

            Eigen::Vector3d b1 = (b - a).normalized();
            Eigen::Vector3d b2 = (c - a).normalized();
            Eigen::Vector3d n = b1.cross(b2).normalized();

            cur_facet_basis.row(0) = b1;
            cur_facet_basis.row(1) = b1.cross(n); //b2


            if (std::abs(n[2]) > .999)
            {
                rot_facet_to_template.setIdentity();
            }
            else
            {
                Eigen::Vector3d c1 = n.cross(b1).normalized();
                
                rot_facet_to_template.col(0) = c1;
                rot_facet_to_template.col(1) = n.cross(c1).normalized();
                rot_facet_to_template.col(2) = n;

            }

            rot_facets_to_template.at(idx) = (rot_facet_to_template.transpose());
            cur_tet_bases.at(idx) = cur_facet_basis;
        }
        appState->R_facet_to_template.at(tetidx) = rot_facets_to_template;
        appState->tet_facet_basis.at(tetidx) = cur_tet_bases;
    }




    /* 

    int nedges = appState->cur_surf->nEdges();
    appState->C_primal.resize(nedges, 2);
    appState->C_sym_2.resize(nedges, 4);
    appState->C_sym_4.resize(nedges, 16); // todo add metric to make in reduced form.  

    // e_projs2.resize(nedges,4); 


    std::vector<Eigen::Matrix2d> rots;// 
    std::vector<Eigen::Matrix4d> rstars;
    std::vector<Eigen::Vector4d> e_projs;



    for (int i = 0; i < nedges; i++)
    {
        Eigen::Vector3d estart = appState->V.row(appState->cur_surf->data().edgeVerts(i, 0));
        Eigen::Vector3d eend = appState->V.row(appState->cur_surf->data().edgeVerts(i, 1));
        Eigen::Vector3d edge_dir = (eend - estart).normalized();
        Eigen::Matrix2d e_to_x;
        e_to_x << edge_dir(0), edge_dir(1), -edge_dir(1), edge_dir(0); // Note this rotates the edge into [1,0]
        // std::cout << e_to_x * edge_dir.head(2) << std::endl<< std::endl; // sanity check.
// std::cout << flatten(edge_dir.head(2) * edge_dir.head(2).transpose()) << std::endl<< std::endl;

        rots.push_back(e_to_x);

        Eigen::Vector4d e_proj = rstar_xcomp_from_r(e_to_x);
        Eigen::VectorXd e_proj_L4_full = rstar_from_r_L4(e_to_x);

        appState->C_primal.row(i) = edge_dir.head(2);
        appState->C_sym_2.row(i) = e_proj;
        appState->C_sym_4.row(i) = e_proj_L4_full;

        // e_projs.push_back(e_proj);
        // e_projs_primal.push_back(edge_dir.head(2));
        // e_projs2.row(i) = e_proj;

    }

    */
//    for(int i = 0; i < appState->R_facet_to_template.size(); i++)
//    {
//         for(int j = 0; j < 4; j++)
//         {
//             std::cout << "R_facet_to_template " << i << " " << j << std::endl;
//             std::cout << appState->R_facet_to_template.at(i).at(j) << std::endl;
//         }
//    }
}

void Mint3DHook::initializeLogFolder() {
    // Retrieve current time
    auto now = std::chrono::system_clock::now();
    auto date = date::floor<date::days>(now);
    auto ymd = date::year_month_day{ date };
    auto time = date::make_time(std::chrono::duration_cast<std::chrono::milliseconds>(now - date));

    // Construct the log folder path using the mesh name and the current date-time
    std::ostringstream folderStream;
    folderStream << "../../results/" << appState->solveType << "/" << appState->meshName << "_"
        << static_cast<unsigned>(ymd.month()) << "_"
        << static_cast<unsigned>(ymd.day()) << "_"
        << time.hours().count() << "_"
        << time.minutes().count();

    // std::cout << "Log folder path: " << folderStream.str() << std::endl;

    // Store the constructed path in AppState
    appState->logFolderPath = folderStream.str();

    // Create the directory using std::filesystem
    std::filesystem::create_directories(appState->logFolderPath);

    std::filesystem::create_directories(appState->logFolderPath + "/outer_iter_ends");
    std::filesystem::create_directories(appState->logFolderPath + "/outer_iter_starts");


    // Inform the user about the log folder creation
    std::cout << "Log folder created at: " << appState->logFolderPath << std::endl;
}


void Mint3DHook::initializeOtherParameters() {
    // ... implementation for initializing other parameters ...
}


void Mint3DHook::initBoundaryConditions() {}
// void Mint3DHook::initBoundaryConditions() {
//     // Assuming boundary faces are identified in AppState
//     Eigen::MatrixXi K;

//     // Eigen::MatrixXi bound_face_idx = appState->bound_face_idx;

//     Eigen::VectorXi boundaryFaces;
//     igl::on_boundary(appState->F, boundaryFaces, K);

//     appState->bound_face_idx = boundaryFaces;

//     // Initialize boundary conditions
//     for (int i = 0; i < boundaryFaces.size(); ++i) {
//         if (boundaryFaces(i) == 1) { // If face is on the boundary
//             Eigen::RowVector3d centroid = (appState->V.row(appState->F(i, 0)) +
//                 appState->V.row(appState->F(i, 1)) +
//                 appState->V.row(appState->F(i, 2))) / 3.0;

//             if (centroid.norm() < 0.45) { // Custom condition for boundary faces
//                 boundaryFaces(i) = -1; // Mark for special handling or exclusion
//             }
//             else {
//                 // Set frame orientation based on the centroid
//                 Eigen::Vector2d frame = Eigen::Vector2d(centroid.y(), -centroid.x()).normalized();
//                 appState->frames.row(i) = frame;
//             }
//         }
//     }

//     appState->frames_orig = appState->frames;

//     // 

//     // 

// }






void Mint3DHook::updateOptimizationParameters() {
    // Example implementation (needs to be adapted to specific needs)
    if (appState->currentIteration % 10 == 0) {
        // Adjust weights or parameters every 10 iterations
        // appState->identityWeight *= 0.9; // Example: gradually decrease the identity weight
    }
    // Add other conditional parameter adjustments as needed
}


void Mint3DHook::checkAndUpdateConvergence(double decrement, double energy) {
    // Example convergence check
    // if (decrement < appState->convergenceThreshold) {
    //     appState->isConverged = true;
    //     std::cout << "Optimization has converged." << std::endl;
    // }
}

void Mint3DHook::finalizeIteration() {
    // Log current state if necessary
    // Check for additional stopping conditions
    // Example:
    if (appState->currentIteration >= appState->maxIterations) {
        std::cout << "Reached maximum number of iterations." << std::endl;
        this->pause(); // Pause the simulation
    }
}


