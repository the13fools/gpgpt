#include "ImGuiWidgets.h"
#include "imgui.h"
// #include "Polyscope.h"
#include "polyscope/polyscope.h"
#include "FieldView.h"

#include "MyConfig.h"

namespace ImGuiWidgets {

    void ShowOptWeights(AppState& appState)
    {
        MyConfig* c = appState.config;
        ImGui::InputDouble("Smoothness Weight", &c->w_smooth);
        ImGui::InputDouble("S Perp Weight", &c->w_s_perp);
        ImGui::InputDouble("Curl Weight", &c->w_curl);
        ImGui::InputDouble("Bound Weight", &c->w_bound);
    }

    // Function to display a file scrubber in ImGui
    void ShowFileScrubber(int& fileIndex, int minIndex, int maxIndex) {
        ImGui::Text("File Scrubber:");

        // Display a slider to select the current file index
        ImGui::SliderInt("File Index", &fileIndex, minIndex, maxIndex);
    }

    // Function to display checkboxes for field views in ImGui
    void ShowFieldViewCheckboxes(AppState& appState) {
        ImGui::Text("Field Views:");

        for (int i = 0; i < Views::Element_COUNT; ++i) {
            Field_View field = static_cast<Field_View>(i);
            const char* fieldName = fieldViewToString(field).c_str();
            bool* active = &appState.fieldViewActive[field];

            ImGui::Checkbox(fieldName, active);
            ImGui::SameLine();
        }
    }

    // Function to display run information in ImGui
    void ShowRunInfo(AppState& appState) {
        ImGui::Text("Run Info:");

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
            std::string fieldName = fieldViewToString(field);
            // std::vector<double> scalarValues(appState.fieldData[field].begin(), appState.fieldData[field].end());

            // Calculate bounds (10th and 90th percentiles)
            float minBound = appState.fieldBounds[field].lower;
            float maxBound = appState.fieldBounds[field].upper;

            // Add the scalar quantity with bounds to Polyscope
            // auto scalarQ = polyscope::getVolumeMesh("my mesh")->addVertexScalarQuantity(fieldName, scalarValues);
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
        ImGui::Text("Field View Scrubber:");
        ImGui::Text("Current Field: %s", fieldViewToString(currentField).c_str());

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

        // Display file scrubber for selecting files
        // ShowFileScrubber(appState.currentFileIndex, 0, appState.fileList.size() - 1);


        ShowOptWeights(appState);

        // Display checkboxes with min and max sliders for field views
        ShowFieldViewCheckboxesWithSliders(appState);

        // Display a field view scrubber
        ShowFieldViewScrubber(appState, appState.current_element);

 ImGui::ShowDemoWindow(); // TODO remove this

        // Display run information
        ShowRunInfo(appState);

        // Display a plot for run_step_times
        // ShowPlot(appState.run_step_times, "Step Times", 0.0f, 100.0f);

        // Display a plot for run_energy
        // ShowPlot(appState.run_energy, "Energy", 0.0f, 100.0f);

        // Additional GUI elements or widgets can be added here

        // End the main ImGui window
        ImGui::End();
    }


}  // namespace ImGuiWidgets
