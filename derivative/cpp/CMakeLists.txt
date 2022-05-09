cmake_minimum_required(VERSION 3.14)
project(numerical_derivative)

set(SRC_DIR src)
set(INCLUDE_DIR include)

set(SRC ${SRC_DIR}/main.cpp ${SRC_DIR}/derivative.cpp)
set(HEADERS ${INCLUDE_DIR}/derivative.hpp)

add_executable(${PROJECT_NAME} ${SRC} ${HEADERS})
target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})