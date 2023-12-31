cmake_minimum_required(VERSION 3.12)
project(04_adfunc_wrapper)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)	  

add_executable(${PROJECT_NAME} ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ../../src)
# add_subdirectory(../libs/TinyAD)

# file(GLOB FRONTENDFILES src/ADWrapper/*.cpp src/ADWrapper/*.h )
# message(STATUS "FRONTENDFILES: ${FRONTENDFILES}")
# target_sources(${PROJECT_NAME} PUBLIC ${FRONTENDFILES})

# Use C++17
target_compile_definitions(${PROJECT_NAME} PUBLIC SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} gpgpt polyscope TinyAD igl::core )

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
