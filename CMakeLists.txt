cmake_minimum_required(VERSION 3.6)
project(Gravitational_potencial)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wall -pedantic -Wextra")

set(SOURCE_FILES main.cpp model.h model.cpp sortByLatitude.h sortByLatitude.cpp localOrthogonalSystem.cpp localOrthogonalSystem.h matrix.cpp matrix.h algorithms.cpp algorithms.h gravAlgorithms.cpp gravAlgorithms.h)
add_executable(Gravitational_potencial ${SOURCE_FILES})