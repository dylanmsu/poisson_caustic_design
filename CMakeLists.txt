cmake_minimum_required(VERSION 3.12)

# Set the project name and C++ standard
project(YourProjectName)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add your source files here
set(SOURCE_FILES 
    main.cpp
    src/polygon_utils.cpp
    src/solver.cpp
    src/utils.cpp
    src/mesh.cpp
    src/bvh.cpp
    # Add more source files as needed
)

# Specify the path to your OpenCV build directory
set(OpenCV_DIR "C:/opencv/build/x64/vc16/lib")  # Adjust the path accordingly

# Find the OpenCV package
find_package(OpenCV REQUIRED PATHS ${OpenCV_DIR})
set(OPENCV_DISABLE_THREAD_SUPPORT ON CACHE BOOL "Disable multi-threading support in OpenCV")

if(NOT OpenCV_FOUND)
    message(FATAL_ERROR "OpenCV not found. Please set OpenCV_DIR to the path containing OpenCVConfig.cmake or add it to CMAKE_PREFIX_PATH.")
endif()

# Create the executable
add_executable(CausticEngineering ${SOURCE_FILES})

# Link with OpenCV libraries
target_link_libraries(CausticEngineering PUBLIC 
    opencv_core
    opencv_imgcodecs
    opencv_highgui
    opencv_imgproc
)

# Optionally, specify the output directory for the executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
