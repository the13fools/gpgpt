#ifndef FILE_PARSER_H
#define FILE_PARSER_H

#include <string>
#include <vector>
#include <optional>
#include <Eigen/Dense>
#include "Serialization.h" // Include the Serialization header

#include "AppState.h"

// Define the FileType enum
enum class FileType {
    BFRA,
    BMOM,
    OBJ,
    FRA,
    MOM,
};

/**
 * The FileParser class is responsible for parsing files of different types
 * within a given directory and retrieving relevant data.
 */
class FileParser {
public:
    // Constructor with the directory to scan
    explicit FileParser(const std::string& directoryPath);

    std::string getFileWithID(const std::string& prefix, const std::string& extension, int fileId);





/*

    // Parse a file with a given ID
    bool parseFileWithID(Eigen::VectorXd& data, FileType fileType, int fileId);

    // Retrieve the largest file's data
    bool parseLargestFile(Eigen::VectorXd& data, FileType fileType);

    void parseToAppState(AppState* appState, int fileId);

    // // Load an .obj file if found
    // bool parseObjFile(Eigen::MatrixXd& V, Eigen::MatrixXi& F);

    // Update the current directory path
    void setDirectoryPath(const std::string& directoryPath);

    */

    std::string objFilePath;
    // int numIDs = 0;

    int minID = 0;
    int maxID = 0;


private:
    std::string directoryPath;
    std::vector<std::string> bfraFiles;
    // std::vector<std::string> bmomFiles;

    // std::string largestBfraFile;
    // std::string largestBmomFile;

    // Helper functions to scan and sort files
    void scanDirectory();
    void findFileBounds();
};

#endif // FILE_PARSER_H
