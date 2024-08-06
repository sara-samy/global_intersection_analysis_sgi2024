cmake_minimum_required(VERSION 3.14)

# Set project name and version
project(gc_project VERSION 1.0)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include FetchContent to download external dependencies
include(FetchContent)

# Fetch IPC-Toolkit
message("\n\n == CMAKE fetching ipc-toolkit from Github\n")
FetchContent_Declare(
    ipc_toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG v1.3.0  # Replace with the desired version tag
)
FetchContent_MakeAvailable(ipc_toolkit)

# Add polyscope
message("\n\n == CMAKE recursively building Polyscope\n")
add_subdirectory("deps/polyscope")

# Find Eigen
message("\n\n == CMAKE looking for Eigen\n")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Add the executable for your project
add_executable(gc_project src/main.cpp)

# Include directories
target_include_directories(gc_project PUBLIC ${ipc_toolkit_SOURCE_DIR}/include)
# Include settings
target_include_directories(gc_project PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libigl/include")

# Link against Polyscope
target_link_libraries(gc_project polyscope)
# Link against Eigen
target_link_libraries(gc_project Eigen3::Eigen)
# Link against IPC-Toolkit
target_link_libraries(gc_project ipc::toolkit)

# Add any other configurations you need
# Example: Set specific compiler options
if(MSVC)
    target_compile_options(gc_project PRIVATE /W4 /WX)
else()
    target_compile_options(gc_project PRIVATE -w)
endif()
