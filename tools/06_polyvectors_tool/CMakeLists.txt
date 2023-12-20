cmake_minimum_required(VERSION 3.12)
project(06_polyvectors)

# Find Eigen
find_package(Eigen3 3.3 REQUIRED)
# find_package(OpenMP REQUIRED)

file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp ../../src/SurfaceFields/*.h ../../src/SurfaceFields/*.cpp ../../src/SurfaceFields/*.hpp)
file(GLOB FRONTENDFILES src/mint2D_gpgpt_hook/*.cpp src/mint2D_gpgpt_hook/*.h src/mint2D_gpgpt_hook/*.hpp )
# file(GLOB ADWRAPPER ../src/04_adfunc_wrapper/ADWrapper/*.cpp ../src/04_adfunc_wrapper/ADWrapper/*.h )

add_executable(${PROJECT_NAME}  ${SRCFILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ../../src ../../lib src/mint2D_gpgpt_hook)
target_sources(${PROJECT_NAME} PUBLIC ${FRONTENDFILES})
# target_sources(${PROJECT_NAME} PUBLIC ${FRONTENDFILES} ${ADWRAPPER})



# add_subdirectory()
# add_subdirectory(../libs/TinyAD)


# Use C++17
target_compile_definitions(${PROJECT_NAME} PUBLIC SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}")


# https://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html 
# here's a good reason to upgrade to C++ 17
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

target_link_libraries(${PROJECT_NAME} gpgpt polyscope TinyAD igl::core args) # OpenMP::OpenMP_CXX)

# add_dependencies(${PROJECT_NAME}  ${FRONTENDFILES})

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17 )
