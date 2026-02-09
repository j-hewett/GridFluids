cmake_minimum_required(VERSION 3.16)

project(GridFluids LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)

qt_standard_project_setup()

qt_add_executable(GridFluids main.cpp
    Fluid.h
    Fluid.cpp
    GridWidget.h
    .gitignore)

target_link_libraries(GridFluids PRIVATE Qt6::Widgets)

set_target_properties(GridFluids PROPERTIES
    WIN32_EXECUTABLE ON
)
