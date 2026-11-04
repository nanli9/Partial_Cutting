cmake_minimum_required(VERSION 3.22)
project(partial_cut LANGUAGES C CXX)

# ------------------------------------------------------------------------------
# Build type (single-config generators only; Visual Studio/Xcode pick per build)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Safety net for CMake >= 4.0, which refuses projects whose cmake_minimum_required
# is below 3.5 (polyscope's bundled deps currently pass, but only just).
if(NOT DEFINED CMAKE_POLICY_VERSION_MINIMUM)
    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
endif()

# Faster MSVC builds of the (many-file) dependency targets.
if(MSVC)
    add_compile_options(/MP)
endif()

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

# Treat dependency headers as system headers (silences their warnings in our TU).
set(DEP_SYSTEM "")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.25)
    set(DEP_SYSTEM SYSTEM)
endif()

# ------------------------------------------------------------------------------
# Pinned dependencies. Exactly one copy of each of Eigen / fcpw is used:
#  * Eigen   : fetched here; fcpw and zombie are pointed at it.
#  * fcpw    : fetched here WITHOUT its submodules (its bundled eigen and slang-rhi
#              are not downloaded).
#  * zombie  : fetched with only the submodules it needs (pcg32, nanoflann, oneTBB);
#              its bundled fcpw and nanobind are not downloaded, and its own
#              CMakeLists.txt is not used (it would add_subdirectory its fcpw copy).
#  * polyscope: fetched with its own submodules (glfw, imgui, glm, ...).

# Eigen 5.0.1 (release tag, commit bc3b39870ecb690a623a3f49149a358b95c5781d)
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        5.0.1
    GIT_SHALLOW    TRUE
    SOURCE_SUBDIR  _no_cmake_   # header-only: do not run Eigen's CMakeLists.txt
)

# fcpw: commit pinned to exactly the one zombie's submodule points at
FetchContent_Declare(fcpw
    GIT_REPOSITORY https://github.com/rohan-sawhney/fcpw.git
    GIT_TAG        e498ca81fe01ae20ea3482485c90379307008de8
    GIT_SUBMODULES ""            # skip bundled eigen + slang-rhi
    ${DEP_SYSTEM}
)

# zombie: no release tags exist; pinned to the current main commit
FetchContent_Declare(zombie
    GIT_REPOSITORY https://github.com/rohan-sawhney/zombie.git
    GIT_TAG        b232b1c99e9407b42fd10c216c31e9e9ad819dd7
    GIT_SUBMODULES deps/pcg32 deps/nanoflann deps/tbb
    SOURCE_SUBDIR  _no_cmake_    # we define the zombie target ourselves below
)

# polyscope v2.6.1 (release tag, commit 0c3dd68b9851417e6b2b976d347adc3250026122)
FetchContent_Declare(polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        v2.6.1
    GIT_SHALLOW    TRUE
    ${DEP_SYSTEM}
)

# ---- dependency options -------------------------------------------------------
# fcpw: CPU-only, no Enoki SIMD, no GPU/Slang, no bindings/demos/tests
set(FCPW_USE_ENOKI             OFF CACHE BOOL "" FORCE)
set(FCPW_ENABLE_GPU_SUPPORT    OFF CACHE BOOL "" FORCE)
set(FCPW_BUILD_BINDINGS        OFF CACHE BOOL "" FORCE)
set(FCPW_BUILD_DEMO            OFF CACHE BOOL "" FORCE)
set(FCPW_BUILD_TESTS           OFF CACHE BOOL "" FORCE)
# oneTBB (zombie's submodule): library only
set(TBB_TEST                   OFF CACHE BOOL "" FORCE)
set(TBB_EXAMPLES               OFF CACHE BOOL "" FORCE)
set(TBB_BENCH                  OFF CACHE BOOL "" FORCE)
set(TBB_INSTALL                OFF CACHE BOOL "" FORCE)
set(TBB_DOC_EXAMPLES           OFF CACHE BOOL "" FORCE)
set(TBB_STRICT                 OFF CACHE BOOL "" FORCE)  # no -Werror / /WX in TBB
set(TBBMALLOC_BUILD            OFF CACHE BOOL "" FORCE)
set(TBBMALLOC_PROXY_BUILD      OFF CACHE BOOL "" FORCE)
set(TCM_BUILD                  OFF CACHE BOOL "" FORCE)  # unsupported on macOS
# polyscope: real OpenGL3/GLFW backend + the headless mock backend
set(POLYSCOPE_BACKEND_OPENGL3_GLFW ON  CACHE BOOL "" FORCE)
set(POLYSCOPE_BACKEND_OPENGL_MOCK  ON  CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(eigen fcpw zombie polyscope)

# ---- Eigen target (single copy, header-only) ------------------------------------
if(NOT TARGET Eigen3::Eigen)
    add_library(eigen_headers INTERFACE)
    target_include_directories(eigen_headers SYSTEM INTERFACE ${eigen_SOURCE_DIR})
    add_library(Eigen3::Eigen ALIAS eigen_headers)
endif()

# fcpw's INTERFACE target lists its (empty, not downloaded) deps/eigen folder as
# an include dir; give it the real Eigen.
target_link_libraries(fcpw INTERFACE Eigen3::Eigen)

# ---- oneTBB from zombie's submodule ---------------------------------------------
# oneTBB only supports being built as a shared library, but GLFW (inside polyscope)
# puts BUILD_SHARED_LIBS=OFF into the cache; override it for TBB only.
set(BUILD_SHARED_LIBS ON)
add_subdirectory(${zombie_SOURCE_DIR}/deps/tbb ${CMAKE_BINARY_DIR}/_deps/tbb-build EXCLUDE_FROM_ALL)
unset(BUILD_SHARED_LIBS)

# ---- zombie target (header-only) --------------------------------------------------
add_library(zombie INTERFACE)
target_include_directories(zombie SYSTEM INTERFACE
    ${zombie_SOURCE_DIR}/include
    ${zombie_SOURCE_DIR}/deps/pcg32
    ${zombie_SOURCE_DIR}/deps/nanoflann/include)
target_link_libraries(zombie INTERFACE fcpw Eigen3::Eigen TBB::tbb)
target_compile_features(zombie INTERFACE cxx_std_17)
if(WIN32)
    # zombie/fcpw headers '#define NOMINMAX' and '_USE_MATH_DEFINES' themselves, but
    # only after <cmath> may already have been pulled in by Eigen, so M_PI can be
    # missing; define them up front (empty value = identical to the in-header
    # #defines, so no redefinition warning).
    target_compile_definitions(zombie INTERFACE "NOMINMAX=" "_USE_MATH_DEFINES=")
endif()
if(MSVC)
    target_compile_options(zombie INTERFACE /bigobj)
endif()

# ------------------------------------------------------------------------------
# Experiment executable
add_executable(partial_cut
    src/main.cpp
    src/cut_problem.h
    src/cut_problem.cpp
    src/validation.h
    src/validation.cpp
    src/walk_cache.h
    src/walk_cache.cpp
    src/validation_reuse.h
    src/validation_reuse.cpp
    src/validation_xi.h
    src/validation_xi.cpp
    src/coupling_gate.h
    src/coupling_gate.cpp)
target_link_libraries(partial_cut PRIVATE zombie polyscope)
set_target_properties(partial_cut PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF)

# Visual Studio: make "Run" launch partial_cut (not ALL_BUILD), from the project folder
# so results.txt lands next to README.md.
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT partial_cut)
set_target_properties(partial_cut PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

if(MSVC)
    target_compile_options(partial_cut PRIVATE /W3 /utf-8 /Zc:__cplusplus)
else()
    target_compile_options(partial_cut PRIVATE -Wall)
endif()

# Put the oneTBB shared library next to the executable on Windows (no rpath there).
if(WIN32)
    add_custom_command(TARGET partial_cut POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:partial_cut> $<TARGET_FILE_DIR:partial_cut>
        COMMAND_EXPAND_LISTS)
endif()
