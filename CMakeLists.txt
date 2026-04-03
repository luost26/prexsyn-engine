cmake_minimum_required(VERSION 3.28)
include(FetchContent)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_third_party")
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_BINARY_DIR}/_third_party")
# list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/third_party/pkgs")

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

if(USE_SYSTEM_BOOST)
    find_package(Boost REQUIRED)
    message(STATUS "Using system Boost: ${Boost_INCLUDE_DIRS}")
else()
    if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_third_party/include/boost")
        execute_process(
            COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/third_party/boost.sh" "${CMAKE_CURRENT_BINARY_DIR}/_third_party"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
            RESULT_VARIABLE BOOST_INSTALL_RESULT)
        if(NOT BOOST_INSTALL_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to install Boost: ${BOOST_INSTALL_RESULT}")
        endif()
    else()
        message(STATUS "Boost already installed. Skipping installation.")
    endif()
    set(Boost_USE_STATIC_LIBS ON)
endif()

if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/_third_party/include/rdkit")
    execute_process(
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/third_party/rdkit.sh" "${CMAKE_CURRENT_BINARY_DIR}/_third_party"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        RESULT_VARIABLE RDKIT_INSTALL_RESULT)
    if(NOT RDKIT_INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to install RDKit: ${RDKIT_INSTALL_RESULT}")
    endif()
else()
    message(STATUS "RDKit already installed. Skipping installation.")
endif()

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

FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_SHALLOW TRUE
    GIT_TAG v3.0.2)
FetchContent_MakeAvailable(pybind11)

set(SPDLOG_BUILD_PIC ON)
FetchContent_Declare(
    spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG v1.17.0)
FetchContent_MakeAvailable(spdlog)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS csrc/*.cpp)
list(FILTER SOURCES EXCLUDE REGEX ".*_test.cpp$")
list(FILTER SOURCES EXCLUDE REGEX "bind.cpp$")
list(FILTER SOURCES EXCLUDE REGEX "main.cpp$")
add_library(prexsyn_obj OBJECT ${SOURCES})
set_target_properties(prexsyn_obj PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(prexsyn_obj PUBLIC ${RDKit_LIBRARIES} OpenMP::OpenMP_CXX spdlog::spdlog)
target_include_directories(prexsyn_obj PUBLIC ${csv_SOURCE_DIR}/single_include ${spdlog_SOURCE_DIR}/include)
add_executable(main csrc/main.cpp)
target_link_libraries(main PRIVATE prexsyn_obj)

find_package(Python COMPONENTS Interpreter Development)
file(GLOB_RECURSE BIND_SOURCES CONFIGURE_DEPENDS csrc/*bind.cpp)
pybind11_add_module(prexsyn_engine ${BIND_SOURCES})
target_link_libraries(prexsyn_engine PRIVATE prexsyn_obj)
if(CMAKE_BUILD_TYPE STREQUAL "Debug" AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(prexsyn_engine PRIVATE -shared-libasan)
endif()
install(TARGETS prexsyn_engine DESTINATION .)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/prexsyn_engine-stubs")
    install(DIRECTORY prexsyn_engine-stubs DESTINATION .)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_SHALLOW TRUE
        GIT_TAG v1.17.0)
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt
        ON
        CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    include(GoogleTest)

    enable_testing()
    file(GLOB_RECURSE TEST_SOURCES csrc/*_test.cpp)
    add_executable(tests ${TEST_SOURCES})
    target_link_libraries(tests PRIVATE prexsyn_obj GTest::gtest_main)
    gtest_discover_tests(tests)
endif()
