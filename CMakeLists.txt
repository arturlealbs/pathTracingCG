cmake_minimum_required(VERSION 3.10)
project(PathTracing)

set(CMAKE_CXX_STANDARD 11)

# Include directories
include_directories(algelin)
include_directories(src)

# Source files
set(ALGELIN_SOURCES
    algelin/OrthographicCamera/OrthoCam.cpp
    algelin/Vector/Vec.cpp
    algelin/Ray/Ray.h
    algelin/Triangle/Triangle.cpp
    algelin/Triangle/Triangle.h
)

set(SRC_SOURCES
    src/Scene.h
    src/Scene.cpp
    src/PathTracing.h
    src/PathTracing.cpp
    src/main.cpp
)

# Executable
add_executable(PathTracing ${ALGELIN_SOURCES} ${SRC_SOURCES})