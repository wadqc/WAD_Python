cmake_minimum_required(VERSION 2.8)

PROJECT(WADAnalysisModule)

# Find ITK.
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)

SUBDIRS(
app
lib
)

SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
#SET(CMAKE_EXE_LINKER_FLAGS -static CACHE STRING "Flags used by the linker" FORCE)
