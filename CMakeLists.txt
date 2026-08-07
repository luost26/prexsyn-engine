cmake_minimum_required(VERSION 3.28)
include(FetchContent)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(prexsyn_engine LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(PREXSYN_DEPS_DIR "${CMAKE_CURRENT_BINARY_DIR}/_third_party" CACHE PATH
    "Directory containing (or receiving) the prexsyn-engine C++ dependencies")
if(NOT IS_ABSOLUTE "${PREXSYN_DEPS_DIR}")
    cmake_path(ABSOLUTE_PATH PREXSYN_DEPS_DIR BASE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
               NORMALIZE)
endif()
option(USE_SYSTEM_BOOST "Use an existing Boost installation" OFF)
option(USE_SYSTEM_RDKIT "Use an existing RDKit installation" OFF)
set(PREXSYN_OPENMP_ROOT "" CACHE PATH "OpenMP runtime prefix (primarily for Apple Clang)")

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
file(MAKE_DIRECTORY "${PREXSYN_DEPS_DIR}")

if(APPLE AND PREXSYN_OPENMP_ROOT)
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "OpenMP compiler flags" FORCE)
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "OpenMP runtime library name" FORCE)
    set(OpenMP_omp_LIBRARY "${PREXSYN_OPENMP_ROOT}/lib/libomp.dylib" CACHE FILEPATH
        "OpenMP runtime library" FORCE)
    set(OpenMP_CXX_INCLUDE_DIR "${PREXSYN_OPENMP_ROOT}/include" CACHE PATH
        "OpenMP include directory" FORCE)
endif()
find_package(OpenMP REQUIRED COMPONENTS CXX)

if(NOT (USE_SYSTEM_BOOST AND USE_SYSTEM_RDKIT))
    set(PREXSYN_DEPS_COMMAND
        "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/third_party/build.py"
        "${PREXSYN_DEPS_DIR}" "--generator=${CMAKE_GENERATOR}")
    if(CMAKE_GENERATOR_PLATFORM)
        list(APPEND PREXSYN_DEPS_COMMAND "--platform=${CMAKE_GENERATOR_PLATFORM}")
    endif()
    if(CMAKE_GENERATOR_TOOLSET)
        list(APPEND PREXSYN_DEPS_COMMAND "--toolset=${CMAKE_GENERATOR_TOOLSET}")
    endif()
    foreach(PREXSYN_CMAKE_PREFIX IN LISTS CMAKE_PREFIX_PATH)
        list(APPEND PREXSYN_DEPS_COMMAND "--cmake-prefix=${PREXSYN_CMAKE_PREFIX}")
    endforeach()
    if(USE_SYSTEM_BOOST)
        list(APPEND PREXSYN_DEPS_COMMAND --skip-boost)
    endif()
    if(USE_SYSTEM_RDKIT)
        list(APPEND PREXSYN_DEPS_COMMAND --skip-rdkit)
    endif()
    execute_process(
        COMMAND ${PREXSYN_DEPS_COMMAND}
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
        COMMAND_ECHO STDOUT
        RESULT_VARIABLE PREXSYN_DEPS_RESULT)
    if(NOT PREXSYN_DEPS_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to install third-party dependencies: ${PREXSYN_DEPS_RESULT}")
    endif()
endif()

list(PREPEND CMAKE_PREFIX_PATH "${PREXSYN_DEPS_DIR}")
if(NOT USE_SYSTEM_BOOST)
    set(Boost_USE_STATIC_LIBS ON)
endif()
find_package(Boost CONFIG REQUIRED)
message(STATUS "Using Boost: ${Boost_INCLUDE_DIRS}")
find_package(RDKit REQUIRED)
message(STATUS "RDKit_INCLUDE_DIRS: ${RDKit_INCLUDE_DIRS}")
set(RDKit_LIBRARIES)
foreach(RDKIT_COMPONENT RDGeneral GraphMol SmilesParse ChemReactions GeneralizedSubstruct DataStructs
                        Fingerprints MolStandardize Descriptors)
    if(TARGET RDKit::${RDKIT_COMPONENT}_static)
        list(APPEND RDKit_LIBRARIES RDKit::${RDKIT_COMPONENT}_static)
    elseif(TARGET RDKit::${RDKIT_COMPONENT})
        list(APPEND RDKit_LIBRARIES RDKit::${RDKIT_COMPONENT})
    else()
        message(FATAL_ERROR "RDKit did not provide the ${RDKIT_COMPONENT} target")
    endif()
endforeach()

FetchContent_Declare(
    csv
    GIT_REPOSITORY https://github.com/vincentlaucsb/csv-parser.git
    GIT_SHALLOW TRUE
    GIT_TAG 2.5.1)
FetchContent_MakeAvailable(csv)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz)
FetchContent_MakeAvailable(json)

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
target_link_libraries(prexsyn_obj PUBLIC ${RDKit_LIBRARIES} OpenMP::OpenMP_CXX spdlog::spdlog
                                         nlohmann_json::nlohmann_json)
target_include_directories(prexsyn_obj PUBLIC ${csv_SOURCE_DIR}/single_include ${spdlog_SOURCE_DIR}/include)
if(MSVC)
    target_compile_options(prexsyn_obj PRIVATE /W4 $<$<CONFIG:Release>:/O2> $<$<CONFIG:Release>:/fp:fast>)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(AMD64|x86_64)$")
        target_compile_options(prexsyn_obj PRIVATE $<$<CONFIG:Release>:/arch:AVX2>)
    endif()
else()
    target_compile_options(prexsyn_obj PRIVATE -Wall -Wextra $<$<CONFIG:Release>:-O3>
                                                $<$<CONFIG:Release>:-ffast-math>)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(AMD64|amd64|x86_64)$")
        target_compile_options(prexsyn_obj PRIVATE $<$<CONFIG:Release>:-mavx2>
                                                    $<$<CONFIG:Release>:-mfma>)
    endif()
endif()
add_executable(main csrc/main.cpp)
target_link_libraries(main PRIVATE prexsyn_obj)

file(GLOB_RECURSE BIND_SOURCES CONFIGURE_DEPENDS csrc/*bind.cpp)
pybind11_add_module(prexsyn_engine ${BIND_SOURCES})
target_link_libraries(prexsyn_engine PRIVATE prexsyn_obj)
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
