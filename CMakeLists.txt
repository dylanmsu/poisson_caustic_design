cmake_minimum_required(VERSION 3.12)

project(CausticEngineering VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)

# Add your source files here
set(SOURCE_FILES 
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    #src/caustic_design.cpp
    #src/polygon_utils.cpp
    #src/solver.cpp
    #src/utils.cpp
    #src/mesh.cpp
    #src/bvh.cpp
)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Dcimg_display=0")

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(CausticEngineering
        MANUAL_FINALIZATION
        ${SOURCE_FILES}
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET CausticEngineering APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(CausticEngineering SHARED
            ${SOURCE_FILES}
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(CausticEngineering
            ${SOURCE_FILES}
        )
    endif()
endif()

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

target_link_libraries(CausticEngineering PRIVATE Qt${QT_VERSION_MAJOR}::Widgets)

# Qt for iOS sets MACOSX_BUNDLE_GUI_IDENTIFIER automatically since Qt 6.1.
# If you are developing for iOS or macOS you should consider setting an
# explicit, fixed bundle identifier manually though.
if(${QT_VERSION} VERSION_LESS 6.1.0)
  set(BUNDLE_ID_OPTION MACOSX_BUNDLE_GUI_IDENTIFIER com.example.CausticEngineering)
endif()
set_target_properties(CausticEngineering PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS CausticEngineering
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(CausticEngineering)
endif()

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
