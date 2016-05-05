################################################################################
# CMake file for pDCX

include_directories(${MPI_CXX_INCLUDE_PATH})
set(LIBS ${MPI_CXX_LIBRARIES})

set(SOURCES dc.cpp)

add_executable(pDCX ${SOURCES})
target_link_libraries(pDCX ${LIBS})

################################################################################
