#include "FileParser.h"
// #include <experimental/filesystem>
#include <filesystem>

#include <filesystem>

#include <algorithm>

#include <igl/readOBJ.h>

#include "Serialization.h" // Include the Serialization header


// namespace fs = std::experimental::filesystem;
namespace fs = std::filesystem;

// https://stackoverflow.com/a/2072890
inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}



FileParser::FileParser(const std::string& directoryPath)
    : directoryPath(directoryPath) {
    scanDirectory();
    findFileBounds();
    std::cout << "FileParser initialized with directory: " << directoryPath << std::endl;
    std::cout << "Found " << bfraFiles.size() << " BFRA files" << std::endl;
}

void FileParser::scanDirectory() {
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            if (ends_with(filename, ".bfra")) {
                bfraFiles.push_back(entry.path().string());
            }
            // } else if (ends_with(filename, ".bmom")) {
            //     bmomFiles.push_back(entry.path().string());
            // } else if (ends_with(filename, ".obj") && objFilePath.empty()) {
            //     objFilePath = entry.path().string(); // Assuming only one .obj file
            // }
        }
    }
    // findLargestIDFile();
}

std::string FileParser::getFileWithID(const std::string& prefix, const std::string& extension, int fileId) {
    std::string targetFile = prefix + std::to_string(fileId + 100000) + extension;
    // std::cout << targetFile << std::endl;
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            if (filename == targetFile) {
                return entry.path().string();
            }
        }
    }
    return ""; // Return empty string if file not found
}

void FileParser::findFileBounds() {
    auto fileIdComparator = [](const std::string& file1, const std::string& file2) {
        int id1 = std::stoi(file1.substr(file1.find_last_of('_') + 1, 6));
        int id2 = std::stoi(file2.substr(file2.find_last_of('_') + 1, 6));
        return id1 > id2;
    };

    if (!bfraFiles.empty()) {
        std::sort(bfraFiles.begin(), bfraFiles.end(), fileIdComparator);
    }

    // The largest ID file will now be at the end of the sorted list
    if (!bfraFiles.empty()) {
        minID = std::stoi(bfraFiles.front().substr(bfraFiles.front().find_last_of('_') + 1, 6));
        maxID = std::stoi(bfraFiles.back().substr(bfraFiles.back().find_last_of('_') + 1, 6));
        std::cout << "minID: " << minID << " maxID: " << maxID << std::endl;
    }



}



/*


void FileParser::scanDirectory() {
    for (const auto& entry : fs::directory_iterator(directoryPath)) {
        if (entry.is_regular_file()) {
            std::string filename = entry.path().filename().string();
            if (ends_with(filename, ".bfra")) {
                bfraFiles.push_back(entry.path().string());
            } else if (ends_with(filename, ".bmom")) {
                bmomFiles.push_back(entry.path().string());
            } else if (ends_with(filename, ".obj") && objFilePath.empty()) {
                objFilePath = entry.path().string(); // Assuming only one .obj file
            }
        }
    }
    // findLargestIDFile();
}



bool FileParser::parseFileWithID(Eigen::VectorXd& data, FileType fileType, int fileId) {
    // Logic to determine the file name based on the ID and file type
    std::string fileName;
    if (fileType == FileType::BFRA) {
        fileName = "primal_" + std::to_string(fileId) + ".bfra";
    } else if (fileType == FileType::BMOM) {
        fileName = "optvars_" + std::to_string(fileId) + ".bmom";
    } else {
        // Handle other file types if necessary
        return false; // Unsupported file type
    }

    std::string filePath = directoryPath + "/" + fileName;

    switch (fileType) {
        case FileType::BFRA:
        case FileType::BMOM:
            return Serialization::deserializeVector(data, filePath);
        case FileType::FRA:
        case FileType::MOM:
            // return readTextFile(filePath, data);
            break;
        case FileType::OBJ:
            // Handle OBJ file parsing
            break;
        default:
            // Handle other file types if necessary
            return false; // Unsupported file type
    }

    return false; // Return false for unsupported file types (e.g., OBJ)
}

bool FileParser::parseLargestFile(Eigen::VectorXd& data, FileType fileType) {
    std::string filePath;
    switch (fileType) {
        case FileType::BFRA:
            filePath = largestBfraFile;
            break;
        case FileType::BMOM:
            filePath = largestBmomFile;
            break;
        default:
            // Handle other file types if necessary
            return false; // Unsupported file type
    }

    switch (fileType) {
        case FileType::BFRA:
        case FileType::BMOM:
            return Serialization::deserializeVector(data, filePath);
        case FileType::FRA:
        case FileType::MOM:
            // return readTextFile(filePath, data);
            break;
        case FileType::OBJ:
            // Handle OBJ file parsing
            break;
        default:
            // Handle other file types if necessary
            return false; // Unsupported file type
    }

    return false; // Return false for unsupported file types (e.g., OBJ)
}

// bool FileParser::parseObjFile(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
//     if (!objFilePath.empty()) {
//         // Implement the parsing logic for the .obj file using libigl
//         if (igl::readOBJ(objFilePath, V, F)) {
//             return true; // Parsing successful
//         } else {
//             return false; // Parsing failed
//         }
//     }
//     return false; // .obj file not found
// }

void FileParser::setDirectoryPath(const std::string& directoryPath) {
    this->directoryPath = directoryPath;
    bfraFiles.clear();
    bmomFiles.clear();
    objFilePath.clear();
    largestBfraFile.clear();
    largestBmomFile.clear();
    scanDirectory();
}

*/

// Other method implementations as necessary...
