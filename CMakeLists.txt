cmake_minimum_required(VERSION 3.12)
cmake_policy(VERSION 3.12)

project(CausticEngineering)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dcimg_display=0")

# Essential include files to build a node addon,
# you should add this line in every CMake.js based project.
include_directories(${CMAKE_JS_INC})

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

# This line will tell CMake that we're building a shared library
# from the above source files
# named after the project's name
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES} ${CMAKE_JS_SRC})
# This line will give our library file a .node extension without any "lib" prefix
set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "" SUFFIX ".node")
# Essential library files to link to a node addon,
# you should add this line in every CMake.js based project.
target_link_libraries(${PROJECT_NAME} ${CMAKE_JS_LIB})
# Include N-API wrappers
execute_process(COMMAND node -p "require('node-addon-api').include"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE NODE_ADDON_API_DIR
        )
string(REPLACE "\n" "" NODE_ADDON_API_DIR ${NODE_ADDON_API_DIR})
string(REPLACE "\"" "" NODE_ADDON_API_DIR ${NODE_ADDON_API_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${NODE_ADDON_API_DIR})
add_definitions(-DNAPI_VERSION=3)

#find_package(PNG REQUIRED)
#find_package(ZLIB REQUIRED)
#find_package(X11 REQUIRED)

# Create the executable
#add_executable(CausticEngineering ${SOURCE_FILES})

# Include directories for PNG and ZLIB
#target_include_directories(CausticEngineering PUBLIC ${PNG_INCLUDE_DIRS})

# Link with libpng and zlib libraries
#target_link_libraries(CausticEngineering PUBLIC ${PNG_LIBRARIES} ${ZLIB_LIBRARIES})

# Link with X11 libraries
#target_link_libraries(CausticEngineering PUBLIC ${X11_LIBRARIES})

# Install the executable
#install(TARGETS CausticEngineering DESTINATION .)

# Include zlib and libpng libraries in the ZIP archive
#install(FILES ${ZLIB_LIBRARIES} ${PNG_LIBRARIES} DESTINATION .)

#include(InstallRequiredSystemLibraries)
#set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
#set(CPACK_PACKAGE_VERSION_MAJOR "${Tutorial_VERSION_MAJOR}")
#set(CPACK_PACKAGE_VERSION_MINOR "${Tutorial_VERSION_MINOR}")
#set(CPACK_SOURCE_GENERATOR "TGZ")
#include(CPack)
