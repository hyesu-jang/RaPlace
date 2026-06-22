cmake_minimum_required(VERSION 3.16)
project(raplace_light CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

find_package(OpenCV REQUIRED COMPONENTS core imgproc imgcodecs)

add_library(raplace raplace.cpp)
target_link_libraries(raplace PUBLIC ${OpenCV_LIBS})
target_include_directories(raplace PUBLIC ${OpenCV_INCLUDE_DIRS}
                           ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(raplace_demo demo.cpp)
target_link_libraries(raplace_demo PRIVATE raplace)
