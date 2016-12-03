cmake_minimum_required(VERSION 2.6)
project(Grondmechanica)
#if(CMAKE_VERSION VERSION_GREATER 3)
    #cmake_policy(SET CMP0054 NEW)
#endif() 
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(/home/mathias/Documents/Grondmechanica/ext/)
include_directories(/home/mathias/Documents/Grondmechanica/src/)
include_directories(/home/mathias/Documents/nanogui/include/)
include_directories(/home/mathias/Documents/nanogui/ext/glad/include)
include_directories(/home/mathias/Documents/nanogui/ext/eigen)
include_directories(/home/mathias/Documents/nanogui/ext/glfw/include)
include_directories(/home/mathias/Documents/nanogui/ext/nanovg/src)
include_directories(/home/mathias/Documents/nanogui/build/)

file(GLOB GRONDMECHANICA_SRC_FILES "./src/*.cpp")
file(GLOB GRONDMECHANICA_HDR_FILES "./src/*.h")

ADD_EXECUTABLE(Grondmechanica ${GRONDMECHANICA_SRC_FILES} ${GRONDMECHANICA_HDR_FILES} ${OPENGL_LIBRARIES}) 

target_link_libraries(Grondmechanica /home/mathias/Documents/Grondmechanica/dep/libnanogui.so)

