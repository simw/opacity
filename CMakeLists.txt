
PROJECT(opacity3)

cmake_minimum_required(VERSION 2.6)
#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

#add definitions, compiler switches, etc.
# ADD_DEFINITIONS(-Wall -O2)
ADD_DEFINITIONS(-Wall -g)

ADD_LIBRARY(libStore2d src/store2d/store.cpp src/store2d/statisticsmc.cpp src/store2d/convergencetable.cpp)
ADD_LIBRARY(libGlvOld src/glv1/qperpdist.cpp src/glv1/zcolldist.cpp src/glv1/radcalcer.cpp src/glv1/qperparraynew.cpp src/glv1/Function.cpp)
ADD_LIBRARY(libRandoms src/randoms/RandDrand.cpp src/randoms/RandSobol.cpp src/randoms/sobol.cpp)
ADD_LIBRARY(libUtils src/parameters.cpp src/constants.cpp src/progressbar.cpp src/timer.cpp)

#list all source files here
SET_SOURCE_FILES_PROPERTIES(src/main.cpp COMPILE_FLAGS -DOPAC1)
ADD_EXECUTABLE(opacity3 src/main.cpp )
TARGET_LINK_LIBRARIES(opacity3 libStore2d libGlvOld libUtils libRandoms)

ADD_EXECUTABLE(qperptest src/main_qperptest.cpp src/glv1/qperparraynew.cpp src/timer.cpp)

 
