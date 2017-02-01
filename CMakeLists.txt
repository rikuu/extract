project(extract)

cmake_minimum_required(VERSION 2.6)

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
#SET (GATB_CORE_EXCLUDE_TOOLS     1)
#SET (GATB_CORE_EXCLUDE_TESTS     1)
#SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# TOOL
################################################################################

# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# we give the headers directories from :
#       - from project source
#       - from GATB-CORE source
#       - from dsk source
include_directories (SYSTEM src  ${gatb-core-includes})

# htslib
SET(HTSLIB_HOME ${PROJECT_SOURCE_DIR}/../htslib/)
LINK_DIRECTORIES(${HTSLIB_HOME})
include_directories(SYSTEM ${HTSLIB_HOME})

include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/thirdparty/)

# we define the files to be compiled
FILE(GLOB EXTRACT_SOURCES src/main.cpp src/io.cpp src/extract.cpp)
FILE(GLOB EXTRACT_TEST_SOURCES src/io.cpp src/io-test.cpp src/extract-test.cpp)

# we define the artifact to be built: the project binary
add_executable("extract" ${EXTRACT_SOURCES})
add_executable("extract-test" ${EXTRACT_TEST_SOURCES})

# we define which libraries to be linked with project binary
target_link_libraries("extract" ${gatb-core-libraries})
target_link_libraries("extract" hts)
target_link_libraries("extract-test" hts)

# we add a compiler flag
SET(CXX_WARNINGS "-Wall -Weverything -Wno-c++98-compat -Wno-vla -Wno-vla-extension -Wno-c++98-compat-pedantic -pedantic")
SET(GCC_COMPILE_MAP_FLAGS "-std=c++0x -g")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_MAP_FLAGS} ${CXX_WARNINGS}")
