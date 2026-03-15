cmake_minimum_required(VERSION 3.28)
include(FetchContent)
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

FetchContent_Declare(
  csv
  GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser.git
  GIT_SHALLOW TRUE
  GIT_TAG 2.5.1)
FetchContent_MakeAvailable(csv)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS csrc/*.cpp)
list(FILTER SOURCES EXCLUDE REGEX ".*_test.cpp$")
list(FILTER SOURCES EXCLUDE REGEX "binding.cpp$")
list(FILTER SOURCES EXCLUDE REGEX "main.cpp$")
add_library(prexsyn_obj OBJECT ${SOURCES})
set_target_properties(prexsyn_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(prexsyn_obj PUBLIC ${RDKit_LIBRARIES} OpenMP::OpenMP_CXX)
target_include_directories(prexsyn_obj PUBLIC ${csv_SOURCE_DIR}/single_include)
target_compile_definitions(
  prexsyn_obj PUBLIC SPDLOG_HEADER_ONLY SPDLOG_FMT_EXTERNAL_HO FMT_HEADER_ONLY)

add_executable(main csrc/main.cpp)
target_link_libraries(main PRIVATE prexsyn_obj)

find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)
file(GLOB_RECURSE BIND_SOURCES CONFIGURE_DEPENDS csrc/*binding.cpp)
pybind11_add_module(prexsyn_engine ${BIND_SOURCES})
target_link_libraries(prexsyn_engine PRIVATE prexsyn_obj)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(prexsyn_engine PRIVATE -shared-libasan)
endif()
install(TARGETS prexsyn_engine DESTINATION .)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  enable_testing()
  find_package(GTest REQUIRED)
  file(GLOB_RECURSE TEST_SOURCES csrc/*_test.cpp)
  add_executable(tests ${TEST_SOURCES})
  target_link_libraries(tests PRIVATE prexsyn_obj GTest::GTest GTest::Main)
  gtest_discover_tests(tests)
endif()
