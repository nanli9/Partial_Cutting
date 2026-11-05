# XFEM on one nodal ring: the partial-cut enrichment of partial_cut_core (the parent project) in a plain cpu XFEM solid
# (backward euler, newton), shown in Polyscope. The library is xfem/, the example examples/nodal_ring/, the self-tests
# tests/. Built as part of the partial_cut solution (the parent's add_subdirectory): Eigen, polyscope and partial_cut_core
# are the parent's targets.

function(xfem_target_options target)
    set_target_properties(${target} PROPERTIES
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        FOLDER XFEM)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W3 /utf-8 /Zc:__cplusplus /bigobj)
        target_compile_definitions(${target} PRIVATE NOMINMAX _USE_MATH_DEFINES)
    else()
        target_compile_options(${target} PRIVATE -Wall)
    endif()
endfunction()

# ---- the library: no dependency on partial_cut, only on Eigen ------------------------
add_library(xfem STATIC
    xfem/types.hpp
    xfem/tet_mesh.hpp
    xfem/tet4.hpp
    xfem/neo_hookean.hpp
    xfem/xfem_cutter.hpp
    xfem/xfem_cutter.cpp
    xfem/xfem_enrichment.hpp
    xfem/xfem_enrichment.cpp
    xfem/ring_export.hpp
    xfem/ring_export.cpp
    xfem/xfem_system.hpp
    xfem/xfem_system.cpp
    xfem/backward_euler.hpp
    xfem/backward_euler.cpp
    xfem/synthetic_ring.hpp
    xfem/synthetic_ring.cpp)
target_include_directories(xfem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(xfem PUBLIC Eigen3::Eigen)
xfem_target_options(xfem)

# ---- the self-tests -------------------------------------------------------------------
add_executable(xfem_tests tests/xfem_tests.cpp)
target_link_libraries(xfem_tests PRIVATE xfem)
xfem_target_options(xfem_tests)

# ---- the example: the ring's enrichment computed directly by partial_cut_core ----------
add_executable(xfem_nodal_ring
    examples/nodal_ring/nodal_ring.cpp
    examples/nodal_ring/live_ring.hpp
    examples/nodal_ring/live_ring.cpp)
target_link_libraries(xfem_nodal_ring PRIVATE xfem partial_cut_core polyscope)
xfem_target_options(xfem_nodal_ring)

# Where the export files are looked for (--data, and the self-tests' real cases): this project's data/ folder.
target_compile_definitions(xfem_nodal_ring PRIVATE XFEM_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
target_compile_definitions(xfem_tests PRIVATE XFEM_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data")
# Where --write-synthetic puts the stand-in rings: next to the tests, not in data/, which holds partial_cut's exports only.
target_compile_definitions(xfem_tests PRIVATE XFEM_SYNTHETIC_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")

# Next to partial_cut.exe (build/<config>/), so the oneTBB library is found there too.
set_target_properties(xfem_nodal_ring xfem_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
if(WIN32)
    add_custom_command(TARGET xfem_nodal_ring POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:xfem_nodal_ring> $<TARGET_FILE_DIR:xfem_nodal_ring>
        COMMAND_EXPAND_LISTS)
endif()
