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

FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG v2.2.1  # Replace with the desired version tag
)
FetchContent_MakeAvailable(polyscope)

# Add the executable for your project
add_executable(${PROJECT_NAME} src/main.cpp)

# Link against Polyscope
target_link_libraries(${PROJECT_NAME} PUBLIC polyscope)
# Link against IPC-Toolkit
target_link_libraries(${PROJECT_NAME} PUBLIC ipc::toolkit)

# Include directories 
target_include_directories(${PROJECT_NAME} PUBLIC 
    ${ipc_toolkit_SOURCE_DIR}/include
    ${polyscope_SOURCE_DIR}/include
)

# Add any other configurations you need
# Example: Set specific compiler options
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -w)
endif()



# Set project name and version
project(xpbd VERSION 1.0)

# Add the executable for your project
add_executable(${PROJECT_NAME} src/simulate.cpp)

# Link against Polyscope
target_link_libraries(${PROJECT_NAME} PUBLIC polyscope)
# Link against IPC-Toolkit
target_link_libraries(${PROJECT_NAME} PUBLIC ipc::toolkit)

# Include directories 
target_include_directories(${PROJECT_NAME} PUBLIC 
    ${ipc_toolkit_SOURCE_DIR}/include
    ${polyscope_SOURCE_DIR}/include
)

# Add any other configurations you need
# Example: Set specific compiler options
if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
else()
    target_compile_options(${PROJECT_NAME} PRIVATE -w)
endif()
