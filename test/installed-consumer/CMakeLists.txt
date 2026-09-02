cmake_minimum_required(VERSION 3.25)

project(nrgljubljana_installed_consumer LANGUAGES CXX)

find_package(nrgljubljana CONFIG REQUIRED)

add_executable(installed_consumer main.cpp)
target_link_libraries(installed_consumer PRIVATE nrgljubljana::nrgljubljana_c)
