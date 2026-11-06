cmake_minimum_required(VERSION 3.25)

project(RoughVolatilityResearchEngine
    VERSION 0.1.0
    DESCRIPTION "Research foundation for rough-volatility derivatives pricing"
    LANGUAGES CXX
)

set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
include(FetchContent)

option(RVRE_ENABLE_SANITIZERS "Enable AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(RVRE_BUILD_EXAMPLES "Build example executables" ON)
option(RVRE_BUILD_BENCHMARKS "Build benchmark executables" OFF)

if(RVRE_ENABLE_SANITIZERS AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    message(FATAL_ERROR "RVRE_ENABLE_SANITIZERS requires a Clang or GNU compiler")
endif()

add_library(rvre_build_options INTERFACE)
target_compile_features(rvre_build_options INTERFACE cxx_std_20)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(rvre_build_options INTERFACE
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
        -Wsign-conversion
        -Wshadow
        -Wnon-virtual-dtor
    )
endif()

if(RVRE_ENABLE_SANITIZERS)
    target_compile_options(rvre_build_options INTERFACE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
    )
    target_link_options(rvre_build_options INTERFACE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
    )
endif()

add_library(rvre
    src/config/experiment_config.cpp
    src/calibration/calibrator.cpp
    src/calibration/market_data.cpp
    src/experiment.cpp
    src/hedging/experiment.cpp
    src/models/heston.cpp
    src/models/rough_bergomi.cpp
    src/pricing/black_scholes.cpp
    src/pricing/monte_carlo.cpp
    src/random/random_stream.cpp
    src/risk/sensitivities.cpp
    src/simulation/brownian_motion.cpp
    src/simulation/time_grid.cpp
    src/simulation/volterra_process.cpp
)
add_library(rvre::rvre ALIAS rvre)
target_include_directories(rvre
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(rvre PUBLIC rvre_build_options)

add_executable(rvre_experiment src/main.cpp)
target_link_libraries(rvre_experiment PRIVATE rvre::rvre)
target_compile_definitions(rvre_experiment
    PRIVATE RVRE_DEFAULT_CONFIG_PATH="${PROJECT_SOURCE_DIR}/config/default.rvre"
)

if(RVRE_BUILD_EXAMPLES)
    add_executable(rvre_rough_bergomi_example examples/rough_bergomi_example.cpp)
    target_link_libraries(rvre_rough_bergomi_example PRIVATE rvre::rvre)
    add_executable(rvre_price examples/pricing_example.cpp)
    target_link_libraries(rvre_price PRIVATE rvre::rvre)
    add_executable(rvre_calibrate examples/calibrate.cpp)
    target_link_libraries(rvre_calibrate PRIVATE rvre::rvre)
    add_executable(rvre_hedge examples/hedging_experiment.cpp)
    target_link_libraries(rvre_hedge PRIVATE rvre::rvre)
    add_custom_target(rvre_examples
        DEPENDS rvre_experiment rvre_rough_bergomi_example rvre_price
                rvre_calibrate rvre_hedge
        COMMENT "Build the RVRE example executables"
    )
endif()

if(RVRE_BUILD_BENCHMARKS)
    add_executable(rvre_volterra_benchmark benchmarks/volterra_benchmark.cpp)
    target_link_libraries(rvre_volterra_benchmark PRIVATE rvre::rvre)
    add_executable(rvre_pricing_benchmark benchmarks/pricing_benchmark.cpp)
    target_link_libraries(rvre_pricing_benchmark PRIVATE rvre::rvre)
    add_executable(rvre_model_comparison benchmarks/model_comparison.cpp)
    target_link_libraries(rvre_model_comparison PRIVATE rvre::rvre)
endif()

if(BUILD_TESTING)
    set(CATCH_BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG 2b60af89e23d28eefc081bc930831ee9d45ea58b
        GIT_PROGRESS FALSE
    )
    FetchContent_MakeAvailable(Catch2)

    add_executable(rvre_tests
        tests/unit/config_tests.cpp
        tests/unit/core_types_tests.cpp
        tests/unit/random_tests.cpp
        tests/unit/brownian_tests.cpp
        tests/unit/volterra_tests.cpp
        tests/unit/rough_bergomi_tests.cpp
        tests/unit/pricing_tests.cpp
        tests/unit/conventional_models_tests.cpp
        tests/unit/calibration_tests.cpp
        tests/unit/risk_hedging_tests.cpp
    )
    target_link_libraries(rvre_tests PRIVATE rvre::rvre Catch2::Catch2WithMain)
    target_compile_definitions(rvre_tests
        PRIVATE RVRE_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
    )
    add_test(NAME rvre_tests COMMAND rvre_tests)
endif()
