cmake_minimum_required(VERSION 3.28)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(prexsyn_engine LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx2 -mfma -ffast-math")

find_package(
  OpenMP
  COMPONENTS CXX
  REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(Boost_USE_STATIC_LIBS ON)

find_package(RDKit REQUIRED)
message(STATUS "RDKit_INCLUDE_DIRS: ${RDKit_INCLUDE_DIRS}")
set(RDKit_LIBRARIES
    RDKit::RDGeneral_static
    RDKit::GraphMol_static
    RDKit::SmilesParse_static
    RDKit::ChemReactions_static
    RDKit::GeneralizedSubstruct_static
    RDKit::DataStructs_static
    RDKit::Fingerprints_static
    RDKit::MolStandardize_static
    RDKit::Descriptors_static)

add_executable(main csrc/main.cpp)
target_link_libraries(main PRIVATE ${RDKit_LIBRARIES} OpenMP::OpenMP_CXX)
