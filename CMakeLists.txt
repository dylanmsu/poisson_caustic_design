cmake_minimum_required(VERSION 3.12)

project(YourProjectName)

# Add your source files here
set(SOURCE_FILES 
    main.cpp
    src/polygon_utils.cpp
    src/solver.cpp
    src/mesh.cpp
    src/bvh.cpp
    # Add more source files as needed
)

# Create the executable
add_executable(CausticEngineering ${SOURCE_FILES})

# Set C++ standard (change as needed)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

# Add any additional dependencies or libraries
target_link_libraries(CausticEngineering PUBLIC ${OpenCV_LIBS})

# Optionally, you can include additional directories
# include_directories(include)

# Optionally, enable additional compiler options
# add_compile_options(-Wall -Wextra)

# Specify the output directory for the executable (change as needed)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})