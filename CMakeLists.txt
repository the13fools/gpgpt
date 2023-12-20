cmake_minimum_required(VERSION 3.9)
project(gpgpt)


# set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "RelWithDebInfo")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

### Compiler options
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 ) # Emit a compile flags file to support completion engines 


if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # set(CMAKE_CXX_COMPILER /usr/bin/clang++-9)

  # message(${which g++})
  # set(CMAKE_CXX_COMPILER /usr/bin/clang++-9 CACHE PATH "" FORCE) # https://stackoverflow.com/a/17275650
  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra -g3")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wunused-local-typedefs -Wsign-compare -Wbuiltin-macro-redefined")

  # SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wunused-local-typedef -Wsign-compare -Wunused-private-field -Wmacro-redefined")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=undefined -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    message("Setting gcc-specific options")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()


  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
  SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-march=native -O3 -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  # set (CMAKE_SHARED_LINKER_FLAGS "-Wl,-ld64")
  set (CMAKE_SHARED_LINKER_FLAGS "-Wl,-ld_classic")

set(OTHER_LDFLAGS "-Wl,-ld_classic")
# set(OTHER_LDFLAGS "-Wl,-weak_reference_mismatches")



  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)

else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


# set(CMAKE_C_COMPILER "/usr/local/Cellar/llvm/5.0.1/bin/clang")
# set(CMAKE_CXX_COMPILER "/usr/local/Cellar/llvm/5.0.1/bin/clang++")

# OPTION (USE_OpenMP "Use OpenMP to enamble <omp.h>" ON)

# Find OpenMP
if(APPLE) # AND USE_OpenMP)
    set(OPENMP_LIBRARIES "/usr/local/Cellar/llvm/5.0.1/lib")
    set(OPENMP_INCLUDES "/usr/local/Cellar/llvm/5.0.1/include")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
      set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
      set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
      set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
      set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
endif()

# Add polyscope
message("\n\n == Init dependencies \n")

# Eigen
# find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# suite-sparse
message("build with suite-sparse")
find_package(SuiteSparse REQUIRED)
include_directories(${SUITESPARSE_INCLUDE_DIRS})


# Libigl
option(LIBIGL_GLFW                "Build target igl::glfw"                OFF)
option(LIBIGL_PNG                 "Build target igl::png"                 OFF)
option(LIBIGL_IMGUI               "Build target igl::imgui"               OFF)
option(LIBIGL_OPENGL              "Build target igl::opengl"              OFF)
option(LIBIGL_PREDICATES          "Build target igl::predicates"          OFF)
option(LIBIGL_COPYLEFT_COMISO     "Build target igl_copyleft::comiso"     OFF)
option(LIBIGL_COPYLEFT_CORE       "Build target igl_copyleft::core"       ON)
option(LIBIGL_COPYLEFT_TETGEN     "Build target igl_copyleft::tetgen"     OFF)
option(LIBIGL_COPYLEFT_CGAL       "Build target igl_copyleft::cgal"       OFF)
option(LIBIGL_RESTRICTED_TRIANGLE "Build target igl_restricted::triangle" OFF)
include(libigl)

# tinyad
include(tinyad)

# polyscope
include(polyscope)

# json
include (json)

# args
include (args)

# spectra
include(spectra)


file(GLOB SRCFILES src/*.cpp src/*.h src/*.hpp)
file(GLOB SURFACEFIELDS src/SurfaceFields/*.h src/SurfaceFields/*.cpp src/SurfaceFields/*.hpp)
file(GLOB FRONTENDFILES src/gpgpt_frontend/*.cpp src/gpgpt_frontend/*.h src/gpgpt_frontend/*.hpp)
file(GLOB ADWRAPPERFILES src/ADWrapper/*.cpp src/ADWrapper/*.h src/ADWrapper/*.hpp)


add_library(${PROJECT_NAME} ${SRCFILES} ${FRONTENDFILES} ${ADWRAPPERFILES} ${SURFACEFIELDS})
target_compile_definitions(${PROJECT_NAME} PUBLIC SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR} ${spectra_SOURCE_DIR}/include/Spectra)

if(MSVC) # AND USE_OpenMP)
  target_link_libraries(${PROJECT_NAME} igl::core polyscope ${SUITESPARSE_LIBRARIES} nlohmann_json::nlohmann_json TinyAD) # stdc++fs
else() # AND USE_OpenMP)
  target_link_libraries(${PROJECT_NAME} igl::core polyscope ${SUITESPARSE_LIBRARIES} nlohmann_json::nlohmann_json TinyAD stdc++fs) # c++experimental
endif()


set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

install(TARGETS ${PROJECT_NAME})

message(${PROJECT_NAME})
message("====== Adding Tools ======")

add_subdirectory(tools/01_basic_test)
add_subdirectory(tools/02_tiny_ad_test)
add_subdirectory(tools/03_disk_bvp_v1)
add_subdirectory(tools/04_adfunc_wrapper)
# add_subdirectory(tools/04_disk_mint2d_v1)
add_subdirectory(tools/05_disk_gpgpt_v1)
add_subdirectory(tools/06_polyvectors_tool)


# add_executable(${PROJECT_NAME} ${SRC_FILES})
# target_compile_definitions(${PROJECT_NAME} PUBLIC SOURCE_PATH="${CMAKE_CURRENT_SOURCE_DIR}")
# target_link_libraries(${PROJECT_NAME} PUBLIC TinyAD polyscope Eigen3::Eigen igl::core )
# # tinyad needs C++17 (not sure if that means this project does too)
# target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
