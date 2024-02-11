cmake_minimum_required(VERSION 3.12)

# Set the project name and C++ standard
project(YourProjectName)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -g)

# Add your source files here
set(SOURCE_FILES 
    main.cpp
    src/polygon_utils.cpp
    src/solver.cpp
    src/utils.cpp
    src/mesh.cpp
    src/bvh.cpp
)

find_package(PNG REQUIRED)

# Create the executable
add_executable(CausticEngineering ${SOURCE_FILES})

# Include the directory containing png.h
target_include_directories(CausticEngineering PUBLIC ${PNG_INCLUDE_DIRS})

# Link with libpng libraries
target_link_libraries(CausticEngineering PUBLIC ${PNG_LIBRARIES})

# Optionally, specify the output directory for the executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
