cmake_minimum_required (VERSION 3.15)
project(hermes C CXX)

set( CMAKE_VERBOSE_MAKEFILE on)
set( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )
set( CMAKE_C_FLAGS "-O3" )
set( CMAKE_CXX_FLAGS "-O3 -ftree-vectorize -std=c++14 -DTHREADS=16" )

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(CMAKE_BUILD_TYPE Debug)

set(HERMES_EXTRA_INCLUDES)
set(HERMES_EXTRA_LIBRARIES)

set (CMAKE_CXX_STANDARD 17)
cmake_policy(SET CMP0048 NEW)

# Set default build-type to release to enable performance improvements
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    #set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")

if(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--as-needed")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--as-needed")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--as-needed")
    message(STATUS "Use --as-needed linker flags!")
endif(CMAKE_COMPILER_IS_GNUCXX AND NOT APPLE)

# ----------------------------------------------------------------------------
# Fix Apple RPATH
# ----------------------------------------------------------------------------
set(CMAKE_MACOSX_RPATH 1)
option(USE_ABSOLUTE_RPATH "Add absolute rpath to all libraries and executables" ON)
if(APPLE OR USE_ABSOLUTE_RPATH)
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(ABSOLUTE_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
    if(NOT IS_ABSOLUTE ${ABSOLUTE_RPATH})
        set(ABSOLUTE_RPATH ${CMAKE_BINARY_DIR}/${ABSOLUTE_RPATH})
    endif()

    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${ABSOLUTE_RPATH}" isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
        message(STATUS "Use absolute RPATH ${ABSOLUTE_RPATH}")
        set(CMAKE_INSTALL_RPATH "${ABSOLUTE_RPATH}")
    endif("${isSystemDir}" STREQUAL "-1")
endif()

# Fix for linker error on mac, Issue: #147
if(APPLE)
    add_definitions(-arch x86_64)
endif(APPLE)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# GSL required dependency
find_package(GSL REQUIRED)
if(GSL_FOUND)
	list(APPEND HERMES_EXTRA_INCLUDES ${GSL_INCLUDE_DIRS})
	list(APPEND HERMES_EXTRA_LIBRARIES ${GSL_LIBRARIES})
endif(GSL_FOUND)

# FFTW3F (optional for turbulent magnetic fields)
find_package(FFTW3F)
if(FFTW3F_FOUND)
    list(APPEND HERMES_EXTRA_INCLUDES ${FFTW3F_INCLUDE_DIR})
    list(APPEND HERMES_EXTRA_LIBRARIES ${FFTW3F_LIBRARIES})
    add_definitions(-DHERMES_HAVE_FFTW3F)
    list(APPEND HERMES_SWIG_DEFINES -DHERMES_HAVE_FFTW3F)
endif(FFTW3F_FOUND)

option(ENABLE_SYS_CFITSIO "System CFITSIO for FITS output" ON)
if(ENABLE_SYS_CFITSIO)
	#set(CFITSIO_EXE_LINKER_FLAGS "-static")
	find_package(CFITSIO REQUIRED)
	if(CFITSIO_FOUND)
        list(APPEND HERMES_EXTRA_INCLUDES ${CFITSIO_INCLUDE_DIR})
        list(APPEND HERMES_EXTRA_LIBRARIES ${CFITSIO_LIBRARIES})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFITSIO_C_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CFITSIO_EXE_LINKER_FLAGS}")
		add_definitions(-DHERMES_HAVE_CFITSIO)
	endif(CFITSIO_FOUND)
endif(ENABLE_SYS_CFITSIO)

find_package(HERMES)
if(HERMES_FOUND)
    list(APPEND HERMES_EXTRA_INCLUDES ${HERMES_INCLUDE_DIR})
    list(APPEND HERMES_EXTRA_LIBRARIES ${HERMES_LIBRARIES})
else() 
    message(FATAL_ERROR "HERMES not found.")	
endif(HERMES_FOUND)

# copy data files in build
file(GLOB HERMES_DATA "data/*fits*")
file(COPY ${HERMES_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# Print out the include and library paths
message(STATUS "HERMES_INCLUDES: ${HERMES_EXTRA_INCLUDES}")
message(STATUS "HERMES_LIBRARIES: ${HERMES_EXTRA_LIBRARIES}")
# Print out definitions
message(STATUS "HERMES_DEFINES: ${HERMES_SWIG_DEFINES}")
message(STATUS "HERMES_CFLAGS: ${CMAKE_C_FLAGS}")
message(STATUS "HERMES_CXXFLAGS: ${CMAKE_CXX_FLAGS}")	
message(STATUS "HERMES_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "HERMES_MODULE_LINKER_FLAGS: ${CMAKE_MODULE_LINKER_FLAGS}")
message(STATUS "HERMES_SHARED_LINKER_FLAGS: ${CMAKE_SHARED_LINKER_FLAGS}")
message(STATUS "HERMES_CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
message(STATUS "HERMES_CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "HERMES_CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "HERMES_CMAKE_INSTALL_RPATH: ${CMAKE_INSTALL_RPATH}")
message(STATUS "HERMES_CMAKE_INSTALL_RPATH_USE_LINK_PATH: ${CMAKE_INSTALL_RPATH_USE_LINK_PATH}")
message(STATUS "HERMES_CMAKE_SKIP_BUILD_RPATH: ${CMAKE_SKIP_BUILD_RPATH}")

# ----------------------------------------------------------------------------
# Playground
# ----------------------------------------------------------------------------
include_directories(include ${HERMES_EXTRA_INCLUDES})

#add_executable(radiomaps src/radiomaps.cpp)
#target_link_libraries(radiomaps ${HERMES_EXTRA_LIBRARIES})

#add_executable(radiospectra src/radiospectra.cpp)
#target_link_libraries(radiospectra ${HERMES_EXTRA_LIBRARIES})

#add_executable(gammamaps src/gammamaps.cpp)
#target_link_libraries(gammamaps ${HERMES_EXTRA_LIBRARIES})

#add_executable(gammaspectra-LE src/gammaspectra_LE.cpp)
#target_link_libraries(gammaspectra-LE ${HERMES_EXTRA_LIBRARIES})

add_executable(gammaspectra-HE src/gammaspectra_HE.cpp)
target_link_libraries(gammaspectra-HE ${HERMES_EXTRA_LIBRARIES})