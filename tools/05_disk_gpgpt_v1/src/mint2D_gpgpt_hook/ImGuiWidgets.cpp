#include "ImGuiWidgets.h"
#include "imgui.h"
// #include "Polyscope.h"
#include "polyscope/polyscope.h"
#include "FieldView.h"

#include "MyConfig.h"

namespace ImGuiWidgets {

    void ShowOptWeights(AppState& appState)
    {
        // auto c = appState.config;
        MyConfig& c = *appState.config;
        
        // ImGui::InputDouble("S Perp Weight", &c.w_s_perp, 0.0f, 0.0f, "%e");
        ImGui::InputDouble("Curl Weight", &c.w_curl, 0.0f, 0.0f, "%e");
        ImGui::InputDouble("Bound Weight", &c.w_bound, 0.0f, 0.0f, "%e");
        ImGui::InputDouble("Smoothness Weight", &c.w_smooth, 0.0f, 0.0f, "%e");
        ImGui::InputDouble("Attenuate Weight", &c.w_attenuate, 0.0f, 0.0f, "%e");
        // ImGui::InputDouble("Bound Weight", &c.w_bound);

        bool* show_frames_as_lines = &appState.show_frames_as_lines;
        bool* show_frames = &appState.show_frames;
        std::string show_frames_checkbox = ("draw vectors as lines##cb");
        ImGui::Checkbox("show frames", show_frames);
        ImGui::SameLine();
        ImGui::Checkbox(show_frames_checkbox.c_str(), show_frames_as_lines);

        // ImGui::InputDouble("Vec Size", &appState.gui_vec_size);

        ImGui::SliderFloat("Vec Size", &appState.gui_vec_size, 1e-4f, 10., "%.16f", ImGuiSliderFlags_Logarithmic);


        // ImGui::InputDouble("L4 alpha", &appState.L4_alpha);

    }

    // Function to display a file scrubber in ImGui
    // void ShowFileScrubber(int& fileIndex, int minIndex, int maxIndex) {
    void ShowFileScrubber(AppState& appState)
    {
        ImGui::Text("File Scrubber:");

        // Display a slider to select the current file index

        int* fileIndex = &appState.currentFileID;
        int minIndex = 0; // appState.minFileIndex;
        int maxIndex = std::max(appState.max_saved_index, 10); // appState.fileList.size() - 1;
        int prevFileIndex = *fileIndex;

        // appState.currentFileIndex, 0, appState.fileList.size() - 1
        if (ImGui::SliderInt("File Index", fileIndex, minIndex, maxIndex) )
        {
            std::cout << "file index changed" << std::endl;
            std::cout << "prev file index " << prevFileIndex << " | current file index " << *fileIndex <<  std::endl;
            appState.shouldReload = true;
        }
    }

    // Function to display checkboxes for field views in ImGui
    void ShowFieldViewCheckboxes(AppState& appState) {
        

       bool* log_state_to_file = &appState.shouldLogData; 
       ImGui::Checkbox("log stuff to file", log_state_to_file);
       if ( ImGui::CollapsingHeader("Select Which Views To Log:") ) 
       {

            for (int i = 0; i < Views::Element_COUNT-1; ++i) {
                Field_View field = static_cast<Field_View>(i);
                // const char* fieldName = fieldViewToString(field).c_str();
                bool* active = &appState.fieldViewActive[field];
                
                
                std::string field_str = fieldViewToFileStub(field);
                std::string checkbox = (field_str + "##cb");


                // std::cout << fieldName << std::endl;
                ImGui::Checkbox(checkbox.c_str(), active);
                // ImGui::SameLine();
            }

       }
    }

    // Function to display run information in ImGui
    void ShowRunInfo(AppState& appState) {
        ImGui::Text("Run Info:");
        ImGui::Text( appState.solveDescription.c_str() );
        ImGui::Text("Status: ");
        ImGui::SameLine();
        ImGui::Text( appState.solveStatus.c_str() );

        ImGui::Text("");

        // ImGui::Text("Current File: %s", appState.fileList[appState.currentFileIndex].c_str());
        // ImGui::Text("Current Element: %s", fieldViewToString(appState.current_element).c_str());
        ImGui::Text("Current Iteration: %d", appState.currentIteration);
        ImGui::Text("Current File ID: %d", appState.currentFileID);

        // ImGui::Text("Optimization Status");
        // ImGui::Text( appState.solveStatus.c_str() );

        // ImGui::Text("Config State:  bound %.0f", (float) appState.config->w_bound);

        // std::cout << "blah" << std::endl;
        //  smooth %f smooth primal %f curl %f", , appState.config->w_smooth, appState.config->w_smooth_vector, appState.config->w_curl);

        ImGui::Text("Config State:");
                // ImGui::Text("Current Step Time: %f", appState.currentStepTime);
        ImGui::Text("Current Energy: %f", appState.os->cur_global_objective_val);
        ImGui::Text("Step Progress: %e | Max Grad Norm: %e | Step Time: %e ms", appState.cur_step_progress, appState.cur_max_gradient_norm, appState.cur_step_time);
        ImGui::Text("smooth primal %.1f bound %.1f curl %.1f smooth %.5f ", (float) appState.config->w_smooth_vector, (float) appState.config->w_bound, (float) appState.config->w_curl, (float) appState.config->w_smooth);
        ImGui::Text("Attenuate weight %e, actual smoothness %e", appState.config->w_attenuate, appState.config->w_attenuate*appState.config->w_smooth);
        ImGui::Text("Relative Residual %e, Absolute Residual %e, Convergence EPS %e", appState.cur_rel_residual, appState.cur_abs_residual, appState.convergenceEpsilon);
            // 
        //     std::cout << "prev obj " << f << " | _prev_step_progress: " << _prev_step_progress << " | newton dec: " << _dec << "  | max_gradient_norm: " << _max_gradient_norm << std::endl;

        //   std::cout << "finshed newton step" << " | solve took " << _prev_step_time << " ms " << std::endl;




        ImGui::Text("TODO");
        ImGui::Text("Fix curl operator");
        // ImGui::Text("refactor optzoo to use stencils");
        // ImGui::Text("reset/load from folder");
        //  ImGui::Text("file scrubber");
        // ImGui::Text("offby1finished issue");
        ImGui::Text("log energy in vector, graph here");
       

        // ImGui::Text("Current Frame: %d", appState.currentFrameIndex);
        // ImGui::Text("Current Moment: %d", appState.currentMomentIndex);

        // Add other run information as needed
    }

    // Function to draw a file loader widget in ImGui
    void DrawFileLoader(AppState& appState) {
        ImGui::Text("File Loader:");

        // ImGui::InputText("Directory", &appState.directoryPath);
        if (ImGui::Button("Load Files")) {
            appState.refreshFileLists();
            // appState.setDefaultFileIndices();
        }
    }

    // Function to add field view scalars to Polyscope
    void AddFieldViewScalarsToPolyscope(AppState& appState) {
        for (int i = 0; i < Views::Element_COUNT; ++i) {
            Field_View field = static_cast<Field_View>(i);
            std::string fieldName = fieldViewToFileStub(field);
            // std::vector<double> scalarValues(appState.fieldData[field].begin(), appState.fieldData[field].end());

            // Calculate bounds (10th and 90th percentiles)
            float minBound = appState.fieldBounds[field].lower;
            float maxBound = appState.fieldBounds[field].upper;

            // Add the scalar quantity with bounds to Polyscope
            // auto scalarQ = polyscope::curls_primal("my mesh")->addVertexScalarQuantity(fieldName, scalarValues);
            // scalarQ->setMapRange({minBound, maxBound});
        }
    }

    // Function to display field view checkboxes with sliders in ImGui
    void ShowFieldViewCheckboxesWithSliders(AppState& appState) {
        ImGui::Text("Field Views with Sliders:");

        for (int i = 0; i < Views::Element_COUNT - 1; ++i) {
            Field_View field = static_cast<Field_View>(i);
            const char* fieldName = fieldViewToFileStub(field).c_str();
            bool* active = &appState.fieldViewActive[field];
            float* minVal = &appState.fieldBounds[field].lower;
            float* maxVal = &appState.fieldBounds[field].upper;

            std::string field_str = fieldViewToFileStub(field);
            // const char* lower = ("l" + field_str).c_str();
            // const char* upper = ("u" + field_str).c_str();
            // const char* checkbox = ("cb" + field_str).c_str();

            // std::string lower = ("l" + field_str).c_str();
            // std::string upper = ("u" + field_str).c_str();
            // std::string checkbox = ("cb" + field_str).c_str();
            
            std::string lower = ("lower ##" + field_str);
            std::string upper = ("upper##" + field_str);
            std::string checkbox = (field_str + "##cb");


            // std::cout << "lower*****" << lower.c_str() << upper.c_str() << checkbox.c_str() << std::endl;
     
            ImGui::PushItemWidth(150);
            ImGui::InputFloat( lower.c_str(), minVal, 0.01f, .10f, "%.3f");
            ImGui::SameLine();
            ImGui::InputFloat( upper.c_str(), maxVal, 0.01f, .10f, "%.3f");
            ImGui::PopItemWidth();
             ImGui::SameLine();
                   ImGui::Checkbox( checkbox.c_str(), active);

            //                    std::cout << fieldName << std::endl;
           
            // ImGui::SliderFloat((std::string(fieldName) + " (Lower Bound)").c_str(), minVal, 0.0f, 1.0f);
            // ImGui::SameLine();
            // ImGui::SliderFloat((std::string(fieldName) + " (Upper Bound)").c_str(), maxVal, 0.0f, 1.0f);
        }
    }

    // Function to display a field view scrubber in ImGui
    void ShowFieldViewScrubber(AppState& appState, Field_View& currentField) {
        // ImGui::Text("Field View Scrubber:");


        

        // Display a combo box to select the current field view
        if (ImGui::BeginCombo("Select Field View", fieldViewToString(currentField).c_str())) {
            for (int i = 0; i < Views::Element_COUNT; ++i) {
                Field_View field = static_cast<Field_View>(i);
                bool is_selected = (currentField == field);
                if (ImGui::Selectable(fieldViewToString(field).c_str(), is_selected)) {
                    currentField = field; // Update the current field view
                }
                if (is_selected) {
                    ImGui::SetItemDefaultFocus(); // Set the default selection
                }
            }
            ImGui::EndCombo();
        }
        ImGui::Text("Field String: %s", fieldViewToFileStub(currentField).c_str());


        bool* active = &appState.fieldViewActive[currentField];
        float* minVal = &appState.fieldBounds[currentField].lower;
        float* maxVal = &appState.fieldBounds[currentField].upper;

        std::string field_str = fieldViewToFileStub(currentField);
            
        std::string percentile_lower = ("lower_bound##percentile" + field_str);
        std::string percentile_upper = ("upper_bound##percentile" + field_str);

        ImGui::PushItemWidth(150);
        ImGui::InputFloat( percentile_lower.c_str(), minVal, 0.01f, .10f, "%.3f");
        ImGui::SameLine();
        ImGui::InputFloat( percentile_upper.c_str(), maxVal, 0.01f, .10f, "%.3f");
        ImGui::PopItemWidth();

        std::string override_lower = ("lower_bound##ovr" + field_str);
        std::string overide_upper = ("upper_bound##ovr" + field_str);

        bool* ovr_active = &appState.override_bounds_active;
        float* ovr_minVal = &appState.override_bounds.lower;
        float* ovr_maxVal = &appState.override_bounds.upper;
                
        std::string checkbox = ("abs override##ovr");


        ImGui::PushItemWidth(150);

        if (*ovr_active)
        {
            ImGui::Text("vvv Absolute Bounds vvv");
        }
        else 
        {
            ImGui::Text("^^^ Percentile Bounds ^^^");
        }

        ImGui::SameLine();
        // std::cout << fieldName << std::endl;
        ImGui::Checkbox(checkbox.c_str(), ovr_active);
        ImGui::InputFloat( override_lower.c_str(), ovr_minVal, 0.001f, .10f, "%.8f");
        ImGui::SameLine();
        ImGui::InputFloat( overide_upper.c_str(), ovr_maxVal, 0.001f, .10f, "%.8f");
        ImGui::PopItemWidth();

        ImGui::PushItemWidth(300);
        ImGui::SliderFloat("override max (log slider)", ovr_maxVal, 1e-16f, 1e-3f, "%.16f", ImGuiSliderFlags_Logarithmic);
        ImGui::PopItemWidth();

    // //  const char* element_names[Field_View::Element_COUNT] = { "Vector Norms", "Delta Norms", "Vector Dirichlet", "Symmetric Dirichlet", "Vector Curl", "Symmetric Curl", "free" };
    //         // const char* current_element_name = (current_element >= 0 && current_element < Field_View::Element_COUNT) ? element_names[current_element] : "Unknown";
    //         const char* current_element_name = fieldViewToString(currentField).c_str();
    //         ImGui::PushItemWidth(300);
    //         ImGui::SliderInt("Shading Mode", (int *) currentField, 0, Views::Field_View::Element_COUNT - 1, current_element_name);
    //         ImGui::PopItemWidth();
        
    }

    // Function to display a plot in ImGui
    void ShowPlot(const char* label, const std::vector<float>& values, float minY, float maxY) {
        ImGui::Text("Plot: %s", label);
        ImGui::PlotLines(label, &values[0], static_cast<int>(values.size()), 0, NULL, minY, maxY, ImVec2(0, 80));
    }

    void ShowMainGUI(AppState& appState) {
        // Begin the main ImGui window
        ImGui::Begin("Main Window");

    
        
        

        // Display checkboxes with min and max sliders for field views
        // ShowFieldViewCheckboxesWithSliders(appState);
    

    // Display file scrubber for selecting files
        ShowFileScrubber(appState);

        // Display a field view scrubber
        ShowFieldViewScrubber(appState, appState.current_element);

        ShowFieldViewCheckboxes(appState);

        ImGui::Begin("Optimization Params");
        ShowOptWeights(appState);
        ImGui::End();

        // Display run information
         ImGui::Begin("Optimization State");
        ShowRunInfo(appState);
        ImGui::End();

         ImGui::ShowDemoWindow(); // TODO remove this

        // Display a plot for run_step_times
        // ShowPlot(appState.run_step_times, "Step Times", 0.0f, 100.0f);

        // Display a plot for run_energy
        // ShowPlot(appState.run_energy, "Energy", 0.0f, 100.0f);

        // Additional GUI elements or widgets can be added here

        // End the main ImGui window
        ImGui::End();
    }


}  // namespace ImGuiWidgets
