cmake_minimum_required (VERSION 3.9)

project (hivision 
  VERSION 0.1.0
  LANGUAGES CXX
)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(CMAKE_POLICY_DEFAULT_CMP0054 NEW)

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_package(PNG REQUIRED)
find_package(hiredis REQUIRED)
find_package(GDAL REQUIRED)
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS"${CMAKE_EXE_LINKER_FLAGS}${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(hivision
  src/server.cpp
  src/vision.cpp
)

include_directories(${PROJECT_SOURCE_DIR}/include/ ${MPI_INCLUDE_PATH} ${PNG_INCLUDE_PATH} ${HIREDIS_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS})

add_definitions(-std=c++11)
ADD_DEFINITIONS(-Wall -O3)

add_executable(hivision_engine ${PROJECT_SOURCE_DIR}/src/vision.cpp)
target_link_libraries(hivision_engine ${GDAL_LIBRARIES} ${HIREDIS_LIBRARIES} ${MPI_LIBRARIES})
install(TARGETS hivision_engine RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})


add_executable(hivision_server ${PROJECT_SOURCE_DIR}/src/server.cpp)
target_link_libraries(hivision_server ${GDAL_LIBRARIES} ${HIREDIS_LIBRARIES} ${Boost_LIBRARIES} ${PNG_LIBRARIES})
install(TARGETS hivision_server RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX})