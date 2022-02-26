cmake_minimum_required(VERSION 3.10)

message("CMAKE_BUILD_TYPE is: '${CMAKE_BUILD_TYPE}'")

set(CMAKE_CXX_FLAGS "-Wall -pedantic --std=c++11 -g")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH cmake_modules)

project(mcrt)

add_library(mcrt_classes photon.cpp detector.cpp layer.cpp slab.cpp 
    phase_functions.cpp)
add_executable(mcrt mcrt.cpp)
target_link_libraries(mcrt mcrt_classes)

enable_testing()
add_executable(mcrt_unit_tests mcrt_unit_tests.cpp)
target_link_libraries(mcrt_unit_tests mcrt_classes)
add_executable(mcrt_integration_tests mcrt_integration_tests.cpp)
target_link_libraries(mcrt_integration_tests mcrt_classes)

# Unit tests
add_test(test_layer_explicit_phase mcrt_unit_tests 1)
add_test(test_make_slab mcrt_unit_tests 2)
add_test(test_make_slab_one_layer mcrt_unit_tests 3)
add_test(test_make_slab_multi_layer mcrt_unit_tests 4)
add_test(test_make_detector_2d mcrt_unit_tests 5)
add_test(test_make_detector_3d mcrt_unit_tests 6)
add_test(test_detect_photon_intensity mcrt_unit_tests 7)
add_test(test_detect_photon_flux mcrt_unit_tests 8)
add_test(test_make_photon_2d mcrt_unit_tests 9)
add_test(test_make_photon_3d mcrt_unit_tests 10)
add_test(test_photon_escp mcrt_unit_tests 11)

# Integration tests
add_test(test_two_layer_angld_pht_mtch_ri mcrt_integration_tests 1)
add_test(test_two_layer_angld_det_mtch_ri_flip mcrt_integration_tests 2)
add_test(test_two_layer_angld_pht_diff_ri mcrt_integration_tests 3)
add_test(test_two_layer_angld_det_diff_ri_flip mcrt_integration_tests 4)
add_test(test_two_layer_angld_pht_and_det_phi_rot0 mcrt_integration_tests 5)
add_test(test_two_layer_angld_pht_and_det_phi_rot180 mcrt_integration_tests 6)
add_test(test_two_layer_angld_det_and_pht_phi_rot180 mcrt_integration_tests 7)
add_test(test_two_layer_thin_slab_intrfc_fluxes_trnsmn mcrt_integration_tests 8)
add_test(test_two_layer_thin_slab_intrfc_edge_fluxes_tos_rad_diff_ri
    mcrt_integration_tests 9)
add_test(test_two_diff_ri_reciprocity mcrt_integration_tests 10)
add_test(test_mltylyr_same_ri mcrt_integration_tests 11)

find_package(OpenMP)
if(OpenMP_CXX_FOUND AND CMAKE_BUILD_TYPE STREQUAL "Release")
    target_link_libraries(mcrt_classes OpenMP::OpenMP_CXX)
    target_link_libraries(mcrt OpenMP::OpenMP_CXX)
    target_link_libraries(mcrt_integration_tests OpenMP::OpenMP_CXX)
endif()

if(CMAKE_COMPILER_IS_GNUCXX AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CodeCoverage)
    append_coverage_compiler_flags()
    setup_target_for_coverage_gcovr_html(NAME cov
        # see https://github.com/bilke/cmake-modules/issues/9
        EXECUTABLE ctest || exit 0
        DEPENDENCIES
        mcrt_unit_tests
        mcrt_integration_tests
        EXCLUDE mcrt.cpp
        )
endif()
