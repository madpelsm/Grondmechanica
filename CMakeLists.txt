cmake_minimum_required(VERSION 2.6)
project(Grondmechanica)
#if(CMAKE_VERSION VERSION_GREATER 3)
    #cmake_policy(SET CMP0054 NEW)
#endif() 
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake/Modules/)

find_package(NANOGUI REQUIRED)

message("nanogui dir  = ${NANOGUI_INCLUDE_DIR}")

IF(WIN32)
include_directories(C:/Users/mathi/Documents/GitHub/nanogui/include/)
include_directories(C:/Users/mathi/Documents/GitHub/nanogui/ext/glad/include)
include_directories(C:/Users/mathi/Documents/GitHub/nanogui/ext/eigen)
include_directories(C:/Users/mathi/Documents/GitHub/nanogui/ext/glfw/include)
include_directories(C:/Users/mathi/Documents/GitHub/nanogui/ext/nanovg/src)
include_directories(C:/Users/mathi/Documents/GitHub/nanogui/build/)
ELSE()
include_directories(/home/mathias/Documents/nanogui/include/)
include_directories(/home/mathias/Documents/nanogui/ext/glad/include)
include_directories(/home/mathias/Documents/nanogui/ext/eigen)
include_directories(/home/mathias/Documents/nanogui/ext/glfw/include)
include_directories(/home/mathias/Documents/nanogui/ext/nanovg/src)
include_directories(/home/mathias/Documents/nanogui/build/)
ENDIF()

include_directories(${PROJECT_SOURCE_DIR}/ext/)
include_directories(${PROJECT_SOURCE_DIR}/src/)

message("source dir = ${PROJECT_SOURCE_DIR}")

file(GLOB GRONDMECHANICA_SRC_FILES "./src/*.cpp")
file(GLOB GRONDMECHANICA_HDR_FILES "./src/*.h")

ADD_EXECUTABLE(Grondmechanica ${GRONDMECHANICA_SRC_FILES} ${GRONDMECHANICA_HDR_FILES} ${OPENGL_LIBRARIES}) 

IF(WIN32)
target_link_libraries(Grondmechanica PROJECT_SOURCE_DIR/dep/nanogui.lib)
ELSE()
target_link_libraries(Grondmechanica PROJECT_SOURCE_DIR/dep/libnanogui.so)
ENDIF()

