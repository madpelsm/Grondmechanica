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


include_directories(${PROJECT_SOURCE_DIR}/ext/nanogui/include/)
include_directories(${PROJECT_SOURCE_DIR}/ext/nanogui/ext/glad/include)
include_directories(${PROJECT_SOURCE_DIR}/ext/nanogui/ext/eigen)
include_directories(${PROJECT_SOURCE_DIR}/ext/nanogui/ext/glfw/include)
include_directories(${PROJECT_SOURCE_DIR}/ext/nanogui/ext/nanovg/src)
include_directories(${PROJECT_SOURCE_DIR}/ext/nanogui/build/)


include_directories(${PROJECT_SOURCE_DIR}/ext/)
include_directories(${PROJECT_SOURCE_DIR}/src/)

message("source dir = ${PROJECT_SOURCE_DIR}")

file(GLOB GRONDMECHANICA_SRC_FILES "./src/*.cpp")
file(GLOB GRONDMECHANICA_HDR_FILES "./src/*.h")

ADD_EXECUTABLE(Grondmechanica ${GRONDMECHANICA_SRC_FILES} ${GRONDMECHANICA_HDR_FILES} ${OPENGL_LIBRARIES}) 

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

IF(WIN32)
    set(NANOGUI_LIB "${PROJECT_SOURCE_DIR}/dep/nanogui.lib")
ELSE()
    if (LINUX)
        set(NANOGUI_LIB "${PROJECT_SOURCE_DIR}/ext/nanogui/build/libnanogui.so")
    ELSE ()
        set(NANOGUI_LIB "${PROJECT_SOURCE_DIR}/ext/nanogui/build/libnanogui.dylib")
    ENDIF()
ENDIF()

target_link_libraries(Grondmechanica ${NANOGUI_LIB})
