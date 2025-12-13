cmake_minimum_required(VERSION 3.15)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(prexsyn_engine LANGUAGES CXX)
find_program(CLANG_PATH clang)
find_program(CLANGXX_PATH clang++)
set(CMAKE_C_COMPILER "${CLANG_PATH}")
set(CMAKE_CXX_COMPILER "${CLANGXX_PATH}")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wextra -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx2 -mfma -ffast-math")

find_package(RDKit REQUIRED)
message(STATUS "RDKit_INCLUDE_DIRS: ${RDKit_INCLUDE_DIRS}")
set(RDKit_LIBRARIES
    RDKit::RDGeneral
    RDKit::GraphMol
    RDKit::SmilesParse
    RDKit::ChemReactions
    RDKit::GeneralizedSubstruct
    RDKit::DataStructs
    RDKit::Fingerprints
    RDKit::MolStandardize
    RDKit::Descriptors
    RDKit::DistGeometry
    RDKit::DistGeomHelpers
    RDKit::ForceField
    RDKit::ForceFieldHelpers
    RDKit::ShapeHelpers
    RDKit::RDGeometryLib
    RDKit::MolAlign
    RDKit::O3AAlign
    RDKit::MolTransforms
    RDKit::ChemTransforms
    RDKit::MolChemicalFeatures)

find_package(Python3 COMPONENTS Interpreter Development NumPy)
message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
message(STATUS "Python3_INCLUDE_DIRS: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")
message(STATUS "Python3_SOABI: ${Python3_SOABI}")

find_package(Boost COMPONENTS system serialization iostreams python numpy
                              REQUIRED)
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message(STATUS "Available Boost_LIBRARIES: ${Boost_LIBRARIES}")

find_package(OpenMP COMPONENTS CXX REQUIRED)
if(OPENMP_FOUND)
  message(STATUS "OpenMP found.")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(LIB_SOURCES
    csrc/chemspace/chemspace.cpp
    csrc/chemspace/indexer.cpp
    csrc/container/building_block_list.cpp
    csrc/container/reaction_list.cpp
    csrc/detokenizer.cpp
    csrc/featurizer/base.cpp
    csrc/featurizer/fingerprint.cpp
    csrc/featurizer/rdkit_descriptors.cpp
    csrc/featurizer/substructures.cpp
    csrc/featurizer/scaffold.cpp
    csrc/featurizer/synthesis.cpp
    csrc/fingerprints.cpp
    csrc/molops.cpp
    csrc/synthesis.cpp)
set(LIB_LINK_LIBRARIES Boost::system Boost::serialization Boost::iostreams
                       ${RDKit_LIBRARIES})
add_library(libprexsyn OBJECT ${LIB_SOURCES})
set_target_properties(libprexsyn PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
target_link_libraries(libprexsyn PRIVATE ${LIB_LINK_LIBRARIES})
target_compile_definitions(
  libprexsyn PRIVATE SPDLOG_HEADER_ONLY SPDLOG_FMT_EXTERNAL_HO FMT_HEADER_ONLY)

function(add_pymod name)
  if(${ARGC} EQUAL 1)
    set(sources csrc/${name}_bind.cpp)
  else()
    set(sources ${ARGN})
  endif()
  python3_add_library(${name} MODULE WITH_SOABI ${sources}
                      $<TARGET_OBJECTS:libprexsyn>)
  target_link_libraries(${name} PRIVATE ${Python3_LIBRARIES} Boost::python
                                        Boost::numpy ${LIB_LINK_LIBRARIES})
endfunction()

add_pymod(building_block_list csrc/container/building_block_list_bind.cpp)
add_pymod(chemspace csrc/chemspace/chemspace_bind.cpp)
add_pymod(detokenizer)
add_pymod(fingerprints)
add_pymod(indexer csrc/chemspace/indexer_bind.cpp)
add_pymod(molops)
add_pymod(pipeline csrc/pipeline/bind.cpp)
add_pymod(reaction_list csrc/container/reaction_list_bind.cpp)
add_pymod(synthesis)
add_pymod(types)
add_pymod(feature csrc/feature/bind.cpp)
add_pymod(featurizer__base csrc/featurizer/base_bind.cpp)
add_pymod(featurizer__fingerprint csrc/featurizer/fingerprint_bind.cpp)
add_pymod(featurizer__rdkit_descriptors
          csrc/featurizer/rdkit_descriptors_bind.cpp)
add_pymod(featurizer__substructures csrc/featurizer/substructures_bind.cpp)
add_pymod(featurizer__scaffold csrc/featurizer/scaffold_bind.cpp)
add_pymod(featurizer__synthesis csrc/featurizer/synthesis_bind.cpp)

if(BUILD_TESTING)
  enable_testing()
  file(COPY resources DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  set(TEST_SOURCES csrc/synthesis_test.cpp csrc/featurizer/product_test.cpp
                   csrc/pipeline/data_buffer_test.cpp)
  find_package(GTest REQUIRED)

  add_executable(unittest ${TEST_SOURCES} $<TARGET_OBJECTS:libprexsyn>)
  target_link_libraries(unittest PRIVATE GTest::gtest_main
                                         ${LIB_LINK_LIBRARIES})
  include(GoogleTest)
  gtest_discover_tests(unittest)
endif()
