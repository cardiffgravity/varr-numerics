# -*- mode: cmake; conding: utf-8; cmake-tab-width: 4; indent-tabs-mode: nil; -*- vim:fenc=utf-8:ft=tcl:et:sw=4:ts=4:sts=4

# -- package info -----------

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(
    varr-numerics
    LANGUAGES C
    VERSION 0.1.0
)

# use GNU standard paths
include(GNUInstallDirs)

# -- dependencies -----------

include(FindGSL)
include_directories(${GSL_INCLUDE_DIRS})
link_libraries(${GSL_LIBRARIES})

# -- build library ----------

# collect headers
set(VARR_HEADERS
    varr/varr.h
    varr/varr_3_over_4.h
    varr/varr_all.h
    varr/varr_atan.h
    varr/varr_cos.h
    varr/varr_exp.h
    varr/varr_extimer.h
    varr/varr_floor_log2.h
    varr/varr_internal.h
    varr/varr_log.h
    varr/varr_phasor.h
    varr/varr_sin.h
    varr/varr_sixth_root.h
)

# list source code
add_library(
    varr SHARED
    src/varr_3_over_4.c
    src/varr_atan.c
    src/varr_cos.c
    src/varr_exp.c
    src/varr_extimer.c
    src/varr_log.c
    src/varr_phasor.c
    src/varr_sin.c
    src/varr_sixth_root.c
)

# set include path
target_include_directories(varr PUBLIC include)

# add compiler flags
target_compile_options(
    varr
    PUBLIC
    -std=c11
    -O3
    -g
    -mtune=native
    -march=native
    -ffast-math
    -ffinite-math-only
    -Wall
    -mavx
    -D_POSIX_C_SOURCE=200112L
)

# add libraries to link
target_link_libraries(
    varr
    m
    ${GSL_LIBRARIES}
)

# set shared library metadata
set_target_properties(varr PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(varr PROPERTIES SOVERSION 1)
set_target_properties(varr PROPERTIES PUBLIC_HEADER "${VARR_HEADERS}")

# define installation paths
install(
    TARGETS varr
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# -- pkgconfig --------------

# install pkgconfig file
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/varr-numerics.pc.in
    ${CMAKE_SOURCE_DIR}/varr-numerics.pc
    @ONLY
)
install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/varr-numerics.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

# -- build tarball ----------

set(CPACK_PACKAGE_MAJOR ${${PROJECT_NAME}_MAJOR_VERSION})
set(CPACK_PACKAGE_MINOR ${${PROJECT_NAME}_MINOR_VERSION})
set(CPACK_PACKAGE_PATCH ${${PROJECT_NAME}_PATCH_VERSION})
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_MAJOR}.${CPACK_PACKAGE_MINOR}.${CPACK_PACKAGE_PATCH}")

set(CPACK_SOURCE_GENERATOR TXZ)
set(CPACK_SOURCE_PACKAGE_FILE_NAME ${PROJECT_NAME}-${${PROJECT_NAME}_VERSION})
set(CPACK_SOURCE_IGNORE_FILES
    "/.*~$/"
    ".*~$"
    "/.git/"
    "/.gitignore"
    "/build/"
    "/commit-specific/"
)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/varr-numerics.spec.in
    ${CMAKE_SOURCE_DIR}/varr-numerics.spec
    @ONLY
)

include(CPack)
