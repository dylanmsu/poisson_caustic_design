cmake_minimum_required(VERSION 3.12)
cmake_policy(VERSION 3.12)

project(CausticEngineering)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Add include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/extern)

# Add your source files here
set(SOURCE_FILES 
    main.cpp
    src/caustic_design.cpp
    src/polygon_utils.cpp
    src/solver.cpp
    src/utils.cpp
    src/mesh.cpp
    src/bvh.cpp
)

find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)

# Create the executable
add_executable(CausticEngineering ${SOURCE_FILES})

# Include directories for PNG and ZLIB
target_include_directories(CausticEngineering PUBLIC ${PNG_INCLUDE_DIRS})

# Link with libpng and zlib libraries
target_link_libraries(CausticEngineering PUBLIC ${PNG_LIBRARIES} ${ZLIB_LIBRARIES})

# Install the executable
install(TARGETS CausticEngineering DESTINATION .)

# Include zlib and libpng libraries in the ZIP archive
install(FILES ${ZLIB_LIBRARIES} ${PNG_LIBRARIES} DESTINATION .)

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${Tutorial_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${Tutorial_VERSION_MINOR}")
set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)
