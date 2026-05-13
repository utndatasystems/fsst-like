cmake_minimum_required(VERSION 3.14)
project(fsst_benchmark CXX)

# Set C++ standard and enable strict mode
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set base compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -march=native -O3 -fno-stack-protector")

# Add include directories
include_directories(fsst src .)

# Optionally define NDEBUG if in release mode
if(NDEBUG)
    add_compile_definitions(NDEBUG)
endif()

# Add Clang-specific flag if needed
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wno-unqualified-std-cast-call)
endif()

# -- OnPair submodule (header-only public API + small static lib) ---------------
set(ONPAIR_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ONPAIR_NATIVE_ARCH ON CACHE BOOL "" FORCE)
set(ONPAIR_ENABLE_LTO ON CACHE BOOL "" FORCE)
add_subdirectory(third-party/onpair_cpp EXCLUDE_FROM_ALL)

# Source files for the executable
set(SOURCES
    main.cpp
    src/BenchmarkDriver.cpp
    src/FsstWrapper.cpp
    src/OnPairWrapper.cpp
    src/Utility.cpp
    src/StateMachine.cpp
    src/MetaStateMachine.cpp
)

# Add the FSST static library path
link_directories(${CMAKE_SOURCE_DIR}/fsst/build)

# Define the executable
add_executable(main ${SOURCES})

# Link against FSST and OnPair
target_link_libraries(main fsst OnPair::onpair)

# Enable Link Time Optimization (LTO)
set_target_properties(main PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
# Optional: Use onpair's helper to ensure flags are fully applied
onpair_optimize_target(main)