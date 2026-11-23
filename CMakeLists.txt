cmake_minimum_required(VERSION 3.16)

project(bed2graph LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(bed2graph_reader
    src/bed_reader.cpp
)
target_include_directories(bed2graph_reader PUBLIC include)

add_executable(bed2graph
    src/main.cpp
)
target_link_libraries(bed2graph PRIVATE bed2graph_reader)

include(CTest)
if(BUILD_TESTING)
    add_executable(bed_reader_tests
        tests/bed_reader_tests.cpp
    )
    target_link_libraries(bed_reader_tests PRIVATE bed2graph_reader)
    add_test(NAME bed_reader_tests COMMAND bed_reader_tests)
endif()
