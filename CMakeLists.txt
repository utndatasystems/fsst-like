cmake_minimum_required(VERSION 3.14)
project(fsst_benchmark CXX)

# Set C++ standard and enable strict mode
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set base compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -march=native -O3 -fno-stack-protector")

# Add include directories
include_directories(
    fsst
    src
    .
    third-party/token-vldb2026/src
    third-party/token-vldb2026/thirdparty/robin_hood
)

# Optionally define NDEBUG if in release mode
if(NDEBUG)
    add_compile_definitions(NDEBUG)
endif()

# Add Clang-specific flag if needed
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wno-unqualified-std-cast-call)
endif()

set(ONPAIR_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(ONPAIR_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(ONPAIR_INSTALL OFF CACHE BOOL "" FORCE)
find_package(Boost 1.81 QUIET)
if(Boost_FOUND AND NOT TARGET Boost::unordered)
    add_library(Boost::unordered INTERFACE IMPORTED)
    target_include_directories(Boost::unordered INTERFACE ${Boost_INCLUDE_DIRS})
endif()
add_subdirectory(third-party/onpair_cpp)

# Source files for the executable
set(SOURCES
    main.cpp
    src/BenchmarkDriver.cpp
    src/decoders/FsstWrapper.cpp
    src/codecs/TokenizerCodec.cpp
    src/codecs/TokenizerDictionary.cpp
    src/decoders/TokenizerDecoder.cpp
    src/SgttCompressorShim.cpp
    src/Utility.cpp
    src/StateMachine.cpp
    src/MetaStateMachine.cpp
    third-party/token-vldb2026/src/compressor/TokenizerCompressor.cpp
)

# Add the FSST static library path
link_directories(${CMAKE_SOURCE_DIR}/fsst/build)

# Define the executable
add_executable(main ${SOURCES})

# Absolute default path to tokenizer dictionary so execution from build/ works.
target_compile_definitions(main PRIVATE
    SGTT_TOKENIZER_RUNTIME_BIN=1
    SGTT_TOKENIZER_BIN_DIR="${CMAKE_SOURCE_DIR}/third-party/token-vldb2026/src/tokens/openai100kpatched"
)

# Link against FSST
target_link_libraries(main fsst OnPair::onpair)

# Enable Link Time Optimization (LTO)
set_target_properties(main PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
