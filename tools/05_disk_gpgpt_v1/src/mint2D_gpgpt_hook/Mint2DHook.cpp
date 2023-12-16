
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
#include <igl/writeOBJ.h>
#include <igl/on_boundary.h>
#include "date.h"
#include "MyConfig.h"
#include "ADWrapper/ADFuncRunner.h"
#include "Surface.h"
#include "UtilsMisc.h"

#include <thread>

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

    // Update the vertex positions and other data in Polyscope
    // polyscope::getSurfaceMesh("c")->updateVertexPositions(appState->V);
    // polyscope::getSurfaceMesh("c")->updateFaceIndices(appState->F);

////
//////  Here we calculate derived quantities. 

    // Check if need to reload, and do this if needed.
    if (appState->shouldReload) {

        std::cout << "do from file reload" << std::endl;

        pause();
        // while(!this->isPaused())
        // {
        //     // wait for pause
        //     std::this_thread::sleep_for(std::chrono::milliseconds(10));
        // }
        if(!this->isPaused())
        {
            std::cout << "failed to pause" << std::endl;
        }

        initSimulation();

    }

/// 
    appState->os->norms_vec = appState->frames.rowwise().norm();
    appState->os->norms_delta = appState->deltas.rowwise().norm();

    ////// The quantities not explicitly set above are set inside of to_passive calls in optzoo.
    // std::cout << opt->get_current_x() << std::endl;
    Eigen::VectorXd tmp = opt->get_current_x();
    std::cout << tmp.rows() << " " << tmp.cols() << std::endl;

    try{
        double cur_obj = opt->eval_func_at(tmp);

        appState->os->cur_global_objective_val = cur_obj;

        std::cout << "cur_obj " << cur_obj << std::endl;
    }
    catch (std::runtime_error& e)
    {
        std::cout << "optimization crashing when evaluating the objective " << e.what() << std::endl;
    }

    


    std::cout << "update render geometry" << std::endl;

    outputData = new OutputState(*appState->os);

    // clear previous passive_vars



    // renderState = new AppState(*appState);

    outputData->frames.clear();

    int frame_rank = appState->primals_layout.rank;
    int vec_dofs = appState->primals_layout.size / frame_rank;
    int frows = appState->frames.rows();

    if (vec_dofs != 2 )
        std::cout << "wrong number of dofs for mint2d" << std::endl;

    for(int vid = 0; vid < frame_rank; vid++)
    {
        Eigen::MatrixXd vec_cur = Eigen::MatrixXd::Zero(appState->frames.rows(), 3);

        vec_cur << appState->frames.block(0, vid * vec_dofs, frows, vec_dofs),  Eigen::MatrixXd::Zero(frows, 1);
        outputData->frames.push_back(vec_cur);
    }

    // outputData->frames.resize(appState->frames.rows(), 3);
    // outputData->frames << appState->frames, Eigen::MatrixXd::Zero(appState->frames.rows(), 1);


    // 
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

        if(appState->max_saved_index < appState->currentFileID)
        {
            appState->max_saved_index = appState->currentFileID;
        }


        // Additional logging for any other fields in AppState as needed
    }
    appState->LogToFile("curr");


    appState->zeroPassiveVars();

 


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

    if (appState->shouldReload && this->isPaused())
    {
        updateRenderGeometry();
    }

    // Depending on the current element view, render different quantities
    // TODO update this to render smoothness and stripe patterns 


    const char* cur_field = fieldViewToFileStub(appState->current_element).c_str();
    
    // std::cout << cur_field << std::endl;

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
            cur_scalar_quantity = outputData->smoothness_L4;
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
        auto cur_scalar_field = polyscope::getSurfaceMesh("c")->addFaceScalarQuantity(cur_field, cur_scalar_quantity);
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

//     if (appState->current_element != Field_View::gui_free)
//     {
//         // cur_scalar_field->setEnabled(true);
//         // cur_scalar_field->setMapRange(appState->fieldViewActive[appState->current_element]);
//         // (cur_field
// // getBoundsPair
//     }

        // Update other visualization properties based on AppState
        // Example: Vector field visualization
        if (appState->showVectorField) {
            //  std::cout << "show vector field" << std::endl;

            //  std::cout << "appState->frames.rows() " << appState->frames.rows() << std::endl;
            //  std::cout << "renderState->frames.rows() " << renderState->frames.rows() << std::endl;
            //  polyscope::getSurfaceMesh("c");
            
            int num_vecs = outputData->frames.size();
            for(int v = 0; v < num_vecs; v++)
            {

                Eigen::MatrixXd cur_vec = outputData->frames[v];


                double color_shift = (v+1.) * 1.0 / num_vecs;

                auto vectorField = polyscope::getSurfaceMesh("c")->addFaceVectorQuantity("Vector Field " + v, cur_vec);
                vectorField->setVectorColor(glm::vec3(color_shift, 0.7, 0.7));
                auto vectorFieldNeg = polyscope::getSurfaceMesh("c")->addFaceVectorQuantity("Vector Field (negative) " + v, (-1.) * cur_vec);
                vectorFieldNeg->setVectorColor(glm::vec3(color_shift, 0.7, 0.7));

                if(appState->show_frames && appState->show_frames_as_lines)
                {
                    vectorField->setEnabled(true);
                    vectorFieldNeg->setEnabled(true);
                    vectorField->setVectorLengthScale(0.01);
                    vectorFieldNeg->setVectorLengthScale(0.01);

                    vectorField->setVectorRadius(0.001);
                    vectorFieldNeg->setVectorRadius(0.001);

                    
                    
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


            


            

            // auto vectorFieldOrig = polyscope::getSurfaceMesh("c")->addFaceVectorQuantity("Vector Field Orig", appState->frames);
        }

        polyscope::requestRedraw();
              




    }

void Mint2DHook::pause() {
    PhysicsHook::pause();
    appState->keepSolving = true;
    appState->solveStatus = "paused";
    std::cout << "paused" << std::endl;
}
    // Pause the simulati




void Mint2DHook::initSimulation() {
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
    Eigen::MatrixXi F; // Temporary storage for faces

    // Check if .bfra and .bmom files exist

    // bool bfraExists = false;
    // bool bmomExists = false;

    // Load mesh and config
    if (create_new_dir) {
       // Load default mesh and set default config
        std::string default_path = std::string(SOURCE_PATH) + "/../shared/" + appState->meshName + ".obj";
        
        // std::string default_path = "/home/josh/Documents/mint_redux/gpgpt/tools/shared/" + cur_mesh_name + ".obj";

        std::cout << default_path << std::endl;
        if (!igl::readOBJ(appState->objFilePath.value_or(default_path), V, F)) {
            std::cerr << "Failed to load mesh from " << appState->objFilePath.value_or(default_path) << std::endl;
            return;
        }
        appState->cur_surf = std::make_unique<Surface>(V, F);

        // Set default configuration
        appState->config = std::make_unique<MyConfig>(); // Assign default values to the config
        initConfigValues();
        // // Serialize default config to a file
        // if (!Serialization::serializeConfig(appState->config, appState->directoryPath + "/config.json")) {
        //     std::cerr << "Failed to save default config to " << appState->directoryPath + "/config.json" << std::endl;
        // }

    } else {

        // fileParser = new FileParser(appState->directoryPath);
        fileParser = std::make_unique<FileParser>(appState->directoryPath);
        appState->max_saved_index = fileParser->maxID;

        // TODO add parsing for this.
        // appState->config = new MyConfig(); // ADD FILE PARSING NOW! 

        appState->objFilePath = fileParser->objFilePath;

        // TODO: Fix this to load from file.
        std::cout << "fileParser->objFilePath " << fileParser->objFilePath << std::endl;
        if (!igl::readOBJ(fileParser->objFilePath, V, F)) {
            std::cerr << "Failed to load mesh from " << fileParser->objFilePath << std::endl;
            return;
        }
        appState->cur_surf = std::make_unique<Surface>(V, F);

        // appState->shouldReload = true;

        


        // load mesh from file 



        // fileParser->
        // bool bfraExists = fileParser.parseLargestFile((appState->frames), FileType::BFRA);
        // bool bmomExists = fileParser.parseLargestFile((appState->deltas), FileType::BMOM);


            // Deserialize configuration from a file
        // if (!Serialization::deserializeConfig(appState->config, appState->directoryPath + "/config.json")) {
        //     std::cerr << "Failed to load config from " << appState->directoryPath + "/config.json" << std::endl;
        //     // Handle error, possibly set default config
        // }
 
    }

    // fieldViewActive = 

    std::cout << "V.rows() " << V.rows() << " F.rows() " << F.rows() << std::endl;
 
 
    // Set mesh data to AppState
    appState->V = V;
    appState->F = F;

    if (appState->shouldReload)
    {
        
        if ( appState->currentFileID == -1 )
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




    // Initialize other parameters and logging folder
    // initializeOtherParameters();

    if (appState->shouldReload == true)
    {
        std::cout << " reload in progress " << std::endl;

        if ( appState->currentFileID == -1 )
        {
            this->resetAppState();
            appState->currentFileID = appState->max_saved_index;
            appState->config = std::make_unique<MyConfig>(); // Assign default values to the config
            initConfigValues();
        }

        
        loadPrimaryData();
        loadGuiState();
        // 
    }

    if (create_new_dir)
    {
        this->resetAppState();
        initializeLogFolder();
        appState->directoryPath = appState->logFolderPath;
    }

    // save OBJ file 
    if (!igl::writeOBJ(appState->directoryPath + "/" + appState->meshName + ".obj", V, F)) {
        std::cerr << "Failed to save mesh to " << appState->logFolderPath << std::endl;
        // return;
    }


    // Register mesh with Polyscope
    polyscope::registerSurfaceMesh("c", appState->V, appState->F);
    polyscope::getSurfaceMesh("c")->setEdgeWidth(0.6);
    polyscope::view::resetCameraToHomeView();
}

bool Mint2DHook::loadPrimaryData() {
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
        } else {
            std::cerr << "File not found for " << info.prefix << " with ID: " << appState->currentFileID << std::endl;
            success = false;
        }
    }

    std::string config_path = fileParser->getFileWithID("config_", ".json", appState->currentFileID);
    std::cout << "config_path " << config_path << std::endl;
    if (!config_path.empty()) {
        if ( !Serialization::deserializeConfig(*appState->config, config_path) ) {
            std::cerr << "Failed to load data for " << "config_" << " from file: " << config_path << std::endl;
            success = false;
        }
    } else {
        std::cerr << "File not found for " << "config_" << " with ID: " << appState->currentFileID << std::endl;
        success = false;
    }
    




    return success;
}

// MyConfig* conf_new = appState->config;
//         if ( !Serialization::deserializeConfig(*conf_new, config_path) ) {
//             std::cerr << "Failed to load data for " << "config_" << " from file: " << config_path << std::endl;
//             success = false;
//         }
//         else 
//         {
//             appState->config = conf_new;
//         }
//         // std::cout << "conf_new->w_curl " << conf_new->w_curl << std::endl;



bool Mint2DHook::loadGuiState() {
    bool success = true;
    for (int i = 0; i < Views::Field_View::gui_free; ++i) {
        Field_View view = static_cast<Field_View>(i);
        std::string fieldStub = Views::fieldViewToFileStub(view) + "_";
        std::string fieldFile = fileParser->getFileWithID(fieldStub, ".bdat", appState->currentFileID);

        if (!fieldFile.empty()) {
            Eigen::VectorXd data;
            if (!Serialization::deserializeVector(data, fieldFile)) {
                std::cerr << "Failed to load data for " << fieldStub << " from file: " << fieldFile << std::endl;
                success = false;
            } else {
                // Update the correct field in OutputState based on 'view'
                // Example: appState.os->norms_vec = data;
            }
        } else {
            std::cerr << "File not found for " << fieldStub << " with ID: " << appState->currentFileID << std::endl;
            success = false;
        }
    }
    return success;
}




bool Mint2DHook::simulateOneStep() {
    appState->solveStatus = "simulate one step";
    int max_iters = appState->maxIterations;
    int cur_iter = appState->currentIteration;
    double convergence_eps = appState->convergenceEpsilon;
    double cur_obj = opt->eval_func_at(opt->get_current_x());
    if (cur_iter < max_iters && cur_obj > convergence_eps && appState->keepSolving)
    {
        cur_iter++;
        std::cout << std::endl << "*** GLOBAL STEP: " << cur_iter << "*** "  << std::endl;
        std::cout << "DOFS in opt" << opt->_cur_x.rows() << std::endl;
        std::cout << "nvars in opt" << opt->get_num_vars() << std::endl; 
        

        opt->take_newton_step( opt->get_current_x() );
        double rel_res_correction = 1. / (cur_obj);
        appState->cur_rel_residual = std::max(opt->_dec * rel_res_correction, 1e-13);
        appState->cur_abs_residual = opt->_dec;
        appState->cur_max_gradient_norm = opt->_max_gradient_norm;
        appState->cur_step_progress = opt->_prev_step_progress;
        appState->cur_step_time = opt->_prev_step_time;
        // std::cout << "cur_obj: " <<  cur_obj  << " conv_eps: " << convergence_eps << "rel_res_correction: " << rel_res_correction << std::endl;

        // Three conditions for convergence:
        // 1. Relative residual is small
        // 2. Absolute residual is small
        // 3. Gradient norm is small or step progress is negative (converged)
        if (appState->cur_rel_residual  < convergence_eps && appState->cur_abs_residual < 1e-3 && (appState->cur_max_gradient_norm < 1e-4 || opt->_prev_step_progress < 1e-10)) 
        {
            std::cout << "**** Converged current step ****" << std::endl;
            std::cout << "Current Objective is " << opt->get_fval_at_x() << std::endl;

            appState->keepSolving = false;

            // In principle should load next config because an optimziation might have a schedule.
            //  cur_iter = max_iters; // break
            //  adhook.reset_params();
        }

        updateAppStateFromOptState();
        appState->solveStatus = "finished one step";

            

     }
        else if (cur_iter == max_iters) 
        {
            // TINYAD_DEBUG_OUT("Final energy: " << func.eval(opt->get_current_x()));
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
    // appState->currentFileID = 0;
    appState->maxIterations = 9999; // Default maximum iterations
    appState->convergenceEpsilon =  1e-12;// 1e-9;
    appState->outerLoopIteration = 0;

    appState->override_bounds.lower = 0;
    appState->override_bounds.upper = 1e-5;
    // appState->shouldReload = false;



    // Resetting mesh data
    appState->frames.setZero(appState->F.rows(), 2);
    appState->deltas.setZero(appState->F.rows(), 4);
    appState->moments.setZero(appState->F.rows(), 0);

    // Resetting derived quantities
    appState->os->norms_vec.setZero(appState->F.rows());
    appState->os->norms_delta.setZero(appState->F.rows());
    appState->os->curls_primal.setZero(appState->F.rows());
    appState->os->curls_sym.setZero(appState->F.rows());
    appState->os->smoothness_primal.setZero(appState->F.rows());
    appState->os->smoothness_L2.setZero(appState->F.rows());
    appState->os->smoothness_L4.setZero(appState->F.rows());
    appState->os->smoothness_L2x2.setZero(appState->F.rows());

    appState->prev_frame_element = Field_View::Element_COUNT;

    // Reinitialize the boundary conditions if needed
    initBoundaryConditions();

    // Initialize symmetric curl operators 
    initCurlOperators();

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



void Mint2DHook::initCurlOperators() 
{
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
        Eigen::Vector3d estart = appState->V.row(appState->cur_surf->data().edgeVerts(i,0));
        Eigen::Vector3d eend = appState->V.row(appState->cur_surf->data().edgeVerts(i,1));
        Eigen::Vector3d edge_dir = (eend - estart).normalized();
        Eigen::Matrix2d e_to_x;
        e_to_x << edge_dir(0),edge_dir(1),-edge_dir(1),edge_dir(0); // Note this rotates the edge into [1,0]
        // std::cout << e_to_x * edge_dir.head(2) << std::endl<< std::endl; // sanity check.
// std::cout << flatten(edge_dir.head(2) * edge_dir.head(2).transpose()) << std::endl<< std::endl;

        rots.push_back(e_to_x);

        Eigen::Vector4d e_proj = rstar_xcomp_from_r(e_to_x);

        appState->C_primal.row(i) = edge_dir.head(2);
        appState->C_sym_2.row(i) = e_proj;

        // e_projs.push_back(e_proj);
        // e_projs_primal.push_back(edge_dir.head(2));
        // e_projs2.row(i) = e_proj;

    }

}

void Mint2DHook::initializeLogFolder() {
    // Retrieve current time
    auto now = std::chrono::system_clock::now();
    auto date = date::floor<date::days>(now);
    auto ymd = date::year_month_day{date};
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
        // appState->identityWeight *= 0.9; // Example: gradually decrease the identity weight
    }
    // Add other conditional parameter adjustments as needed
}


void Mint2DHook::checkAndUpdateConvergence(double decrement, double energy) {
    // Example convergence check
    // if (decrement < appState->convergenceThreshold) {
    //     appState->isConverged = true;
    //     std::cout << "Optimization has converged." << std::endl;
    // }
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