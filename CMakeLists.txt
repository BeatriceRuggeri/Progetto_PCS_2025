cmake_minimum_required(VERSION 3.0.0)

# ProgettoPCS is the project name
project(ProgettoPCS LANGUAGES C CXX VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND ProgettoPCS_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND ProgettoPCS_LINKED_LIBRARIES PRIVATE Threads::Threads)

################################################################################

## GTest

find_package(GTest REQUIRED)

list(APPEND ProgettoPCS_LINKED_LIBRARIES PRIVATE GTest::gtest)
list(APPEND ProgettoPCS_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
list(APPEND ProgettoPCS_LINKED_LIBRARIES PRIVATE GTest::gmock)
list(APPEND ProgettoPCS_LINKED_LIBRARIES PRIVATE GTest::gmock_main)

################################################################################

# message(${CMAKE_SOURCE_DIR})
# message(${CMAKE_CURRENT_SOURCE_DIR})
# message(${CMAKE_BINARY_DIR})
# message(${CMAKE_CURRENT_BINARY_DIR})
# message(${PROJECT_NAME})
# message(${CMAKE_PROJECT_NAME})

###################### Add subdirectory ##################################

add_subdirectory(src)
add_subdirectory(src_test)
add_subdirectory(ExportParaview)

list(APPEND ProgettoPCS_sources ${src_sources} ${src_paraview_sources})
list(APPEND ProgettoPCS_headers ${src_headers} ${src_paraview_headers})
list(APPEND ProgettoPCS_includes ${src_includes} ${src_paraview_includes})

list(APPEND ProgettoPCS_TEST_sources ${${PROJECT_NAME}_sources} ${src_test_sources})
list(APPEND ProgettoPCS_TEST_headers ${${PROJECT_NAME}_headers} ${src_test_headers})
list(APPEND ProgettoPCS_TEST_includes ${${PROJECT_NAME}_includes} ${src_test_includes})

############################## Add executable ##################################################

add_executable(${PROJECT_NAME} main.cpp
      ${${PROJECT_NAME}_sources}
      ${${PROJECT_NAME}_headers})

target_include_directories(${PROJECT_NAME} PUBLIC ${${CMAKE_PROJECT_NAME}_includes})
target_link_libraries(${PROJECT_NAME} PUBLIC ${${PROJECT_NAME}_LINKED_LIBRARIES})
target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

########################################### executable testing ###############################

add_executable(${PROJECT_NAME}_TEST main_test.cpp
      ${${PROJECT_NAME}_TEST_sources}
      ${${PROJECT_NAME}_TEST_headers}
      ${${PROJECT_NAME}_sources}
      ${${PROJECT_NAME}_headers})

target_include_directories(${PROJECT_NAME}_TEST PUBLIC ${${CMAKE_PROJECT_NAME}_TEST_includes})
target_link_libraries(${PROJECT_NAME}_TEST PUBLIC ${${PROJECT_NAME}_LINKED_LIBRARIES})
target_compile_options(${PROJECT_NAME}_TEST PUBLIC -fPIC)
