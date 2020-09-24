project(KondoHeisenberg)
cmake_minimum_required(VERSION 3.3)
set(CMAKE_CXX_STANDARD 17)

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

set(SOURCES
    main.cpp
    Hamiltonian.cpp
)
set(HEADERS
    Hamiltonian.h
)


add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC include)

target_link_libraries(${PROJECT_NAME} ${GTEST_LIBRARIES} -larmadillo -std=c++17)
