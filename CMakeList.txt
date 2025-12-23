cmake_minimum_required(VERSION 3.15)

# ==============================================================================
# PROJECT METADATA
# ==============================================================================
project(AureliaFinslerSim
    VERSION 1.0.0
    DESCRIPTION "Non-Equilibrium Finslerian Field Theory Simulation Engine"
    LANGUAGES CXX
)

# ==============================================================================
# C++ STANDARD CONFIGURATION (C++17 REQUIRED)
# ==============================================================================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Enforce strict ISO C++

# ==============================================================================
# COMPILER OPTIMIZATIONS (CRITICAL FOR PHYSICS)
# ==============================================================================
# Default to Release mode if not specified (we need speed!)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # -O3           : Maximize optimization
    # -march=native : Use CPU-specific instructions (AVX/AVX2) for matrix math
    # -Wall -Wextra : High warning level (clean code requirement)
    # -funroll-loops: Optimize iteration loops for WLC integration
    add_compile_options(-O3 -march=native -Wall -Wextra -funroll-loops)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # /O2 : Maximize speed
    # /W4 : High warning level
    add_compile_options(/O2 /W4)
endif()

# ==============================================================================
# SOURCE MANAGEMENT
# ==============================================================================
# Locate all header files recursively.
# This ensures that your complicated folder structure (Geometry, ActiveMatter, etc.)
# is recognized by IDEs (VS Code, CLion, etc.) and included in the project tree.
file(GLOB_RECURSE HEADER_FILES
    "config/*.h"
    "src/*.h"
)

# ==============================================================================
# EXECUTABLE TARGET
# ==============================================================================
add_executable(AureliaSim
    main.cpp
    ${HEADER_FILES} # Listed here for IDE visibility, not strictly for compilation
)

# ==============================================================================
# INCLUDE DIRECTORIES
# ==============================================================================
# Allows writing #include "src/Geometry/..." instead of relative paths
target_include_directories(AureliaSim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# ==============================================================================
# POST-BUILD INFO
# ==============================================================================
add_custom_command(TARGET AureliaSim POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "======================================================="
    COMMAND ${CMAKE_COMMAND} -E echo "  AURELIA SIMULATION ENGINE BUILT SUCCESSFULLY"
    COMMAND ${CMAKE_COMMAND} -E echo "  Run with: ./AureliaSim"
    COMMAND ${CMAKE_COMMAND} -E echo "======================================================="
)