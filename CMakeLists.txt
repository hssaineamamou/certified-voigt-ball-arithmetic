
cmake_minimum_required(VERSION 3.16)

project(
    CertifiedVoigt
    VERSION 1.0.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find FLINT through pkg-config when available.
find_package(PkgConfig QUIET)

if(PkgConfig_FOUND)
    pkg_check_modules(FLINT QUIET IMPORTED_TARGET flint)
endif()

add_executable(
    voigt
    voigt.cpp
)

if(TARGET PkgConfig::FLINT)

    target_link_libraries(
        voigt
        PRIVATE
        PkgConfig::FLINT
    )

else()

    find_path(
        FLINT_INCLUDE_DIR
        NAMES flint/flint.h flint/arb.h flint/acb.h
    )

    find_library(
        FLINT_LIBRARY
        NAMES flint
    )

    find_library(
        GMP_LIBRARY
        NAMES gmp
    )

    find_library(
        MPFR_LIBRARY
        NAMES mpfr
    )

    if(NOT FLINT_INCLUDE_DIR OR NOT FLINT_LIBRARY)
        message(FATAL_ERROR
            "FLINT library or headers not found."
        )
    endif()

    target_include_directories(
        voigt
        PRIVATE
        ${FLINT_INCLUDE_DIR}
    )

    target_link_libraries(
        voigt
        PRIVATE
        ${FLINT_LIBRARY}
    )

    if(GMP_LIBRARY)
        target_link_libraries(voigt PRIVATE ${GMP_LIBRARY})
    endif()

    if(MPFR_LIBRARY)
        target_link_libraries(voigt PRIVATE ${MPFR_LIBRARY})
    endif()

endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")

    target_compile_options(
        voigt
        PRIVATE
        -O2
    )

endif()

message(STATUS "Voigt function project configured.")
