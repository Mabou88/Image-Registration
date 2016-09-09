cmake_minimum_required(VERSION 3.2)

project(Image-Registration)
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})
add_executable(Image-Registration main_registration.cpp)
target_link_libraries(Image-Registration ${ITK_LIBRARIES})