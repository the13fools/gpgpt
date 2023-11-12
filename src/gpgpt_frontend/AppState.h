#ifndef APPSTATE_H
#define APPSTATE_H

#include <string>
#include <vector>
#include <optional>
#include <unordered_map>
#include <limits>

#include "FieldView.h"

// // Enum for identifying field view quantities
// enum Field_View {
//     vec_norms,
//     delta_norms,
//     vec_dirch,
//     moment_dirch,
//     primal_curl_residual,
//     sym_curl_residual,
//     gui_free,
//     Element_COUNT // This should always be last
// };

// Struct to hold bounds for each field view quantity
struct FieldBounds {
    float upper = std::numeric_limits<double>::max();
    float lower = std::numeric_limits<double>::lowest();
};

// AppState holds the state of the application
class AppState {
public:
    std::string directoryPath;
    std::vector<std::string> bfraFiles;
    std::vector<std::string> bmomFiles;
    std::optional<std::string> objFilePath;
    int currentFileID = 0;
    std::unordered_map<Field_View, FieldBounds> fieldBounds;
    bool fieldViewActive [];

    // Constructor
    AppState();

    // Method to refresh the file lists
    void refreshFileLists();

    // Define functions expected by GUIContext
    void selectFile(const std::string& filename);
    void refreshData();
    void updateSliderValue(int value);
    void serializeData();
    // More methods as needed
};

#endif // APPSTATE_H
