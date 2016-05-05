################################################################################
# CMake file for pDCX project

project(pDCX CXX)

cmake_minimum_required(VERSION 2.6.2)

enable_testing()

# options


# disallow in-source builds

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# enable warnings

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# use MPI library
find_package(MPI REQUIRED)

# descend into source

add_subdirectory(src)

################################################################################
