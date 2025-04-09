cmake_minimum_required(VERSION 3.16)
project(fsst-like)


# Compiler flags
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include directories
include_directories(. ./src/)

# Files
file(GLOB_RECURSE SOURCE_FILES RELATIVE ${CMAKE_SOURCE_DIR} "src/*/*.[tch]pp") # Source files

# Add entry points (i.e., main functions)
set(MAIN_FILE main.cpp)

add_executable(fsst-like ${MAIN_FILE} ${SOURCE_FILES})
