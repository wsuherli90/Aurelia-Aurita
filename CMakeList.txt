cmake_minimum_required(VERSION 3.15)

# ==============================================================================
# PROJECT METADATA
# ==============================================================================
project(AureliaFinslerSim
    VERSION 1.1.0
    DESCRIPTION "HPC Non-Equilibrium Finslerian Field Theory Simulation Engine"
    LANGUAGES CXX
)

# ==============================================================================
# C++ STANDARD (C++17 Required for std::array & constexpr)
# ==============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ==============================================================================
# HPC COMPILER FLAGS (The "Speed" Sauce)
# ==============================================================================
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # -O3           : Aggressive optimization
    # -march=native : Use all instructions available on THIS cpu (AVX/SSE)
    # -ffast-math   : Allow math optimizations that might violate strict IEEE compliance (faster sqrt, etc)
    add_compile_options(-O3 -march=native -ffast-math -Wall -Wextra -funroll-loops)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # /O2 : Maximize Speed
    # /arch:AVX2 : Enable vectorization
    add_compile_options(/O2 /arch:AVX2 /W4)
endif()

# ==============================================================================
# SOURCE MANAGEMENT
# ==============================================================================
file(GLOB_RECURSE HEADER_FILES
    "config/*.h"
    "src/*.h"
)

# ==============================================================================
# TARGET DEFINITION
# ==============================================================================
add_executable(AureliaSim
    main.cpp
    ${HEADER_FILES} 
)

# Include directories
target_include_directories(AureliaSim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

message(STATUS "HPC Configuration Ready. Run in Release mode for max performance.")