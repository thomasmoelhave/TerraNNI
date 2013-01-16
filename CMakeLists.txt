cmake_minimum_required(VERSION 2.6)
project("TerraNNI-project")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")
include_directories(${CMAKE_SOURCE_DIR})

###############################################
# Setup Dependencies 
###############################################

#####################
# Find TPIE 
option(USE_EXTERNAL_TPIE "Do not use the TPIE that was pulled using git submodules" OFF)
mark_as_advanced(USE_EXTERNAL_TPIE)

if (NOT USE_EXTERNAL_TPIE)
	#
	### We use the TPIE placed in the tpie subdirectory using git submodules
	#
	set(TPIE_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/tpie/)
	message(STATUS "Using TPIE from git submodules at ${CMAKE_SOURCE_DIR}")
	set(TPIE_USE_EXCEPTIONS ON FORCE)
	add_subdirectory(tpie)
	include_directories(${CMAKE_BINARY_DIR}/tpie)
else(NOT USE_EXTERNAL_TPIE)
	find_package(TPIE REQUIRED)
endif(NOT USE_EXTERNAL_TPIE)

include_directories(${TPIE_INCLUDE_DIRS})
link_directories(${TPIE_LIBRARY_DIRS})

#####################
# Find CUDA
find_package(CUDA REQUIRED)
if(CUDA_FOUND)
	set(CUDA_COMPILER_BINDIR "/usr/bin" CACHE PATH "Path for g++ binary for nvcc")
	SET(CUDA_NVCC_FLAGS "-arch=compute_20" "--compiler-bindir=${CUDA_COMPILER_BINDIR}")
	include_directories("${CUDA_SDK_ROOT_DIR}/common/inc/")
	CUDA_ADD_LIBRARY(cuda_nni "cuda_nni.cu")
	set(NNI_LIBS cuda_nni)
endif(CUDA_FOUND)

#####################
# Find Boost 
set(Boost_ADDITIONAL_VERSIONS "1.40.0" "1.40" "1.41" "1.41.0" "1.42" "1.42.0" "1.43" "1.43.0" "1.44" "1.44.0")
set(Boost_USE_STATIC_LIBS    ON)
find_package(Boost COMPONENTS thread system filesystem program_options REQUIRED)
link_directories ( ${Boost_LIBRARY_DIRS} )
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(-DBOOST_FILESYSTEM_VERSION=2)

#####################
# libLAS is required 
find_package(libLAS REQUIRED)
include_directories(${LIBLAS_INCLUDE_DIR})


###############################################
# Setup Libraries 
###############################################
IF(APPLE)
	include_directories(/usr/X11R6/include/)
	link_directories(/usr/X11R6/lib)
ENDIF(APPLE)
SET(EXTRA_LIBS GL X11 GLU glut GLEW)

add_library(nni "nni.cpp")
target_link_libraries(nni  ${EXTRA_LIBS} tpie ${NNI_LIBS})

add_executable(terranni "terranni.cpp")

add_executable(terranni-2d "terranni.cpp")
set_target_properties(terranni-2d PROPERTIES COMPILE_FLAGS -DTERRANNI_2DONLY)

target_link_libraries(terranni  nni tpie ${Boost_LIBRARIES} ${EXTRA_LIBS} ${LIBLAS_LIBRARY} ${TPIE_LIBRARY})
target_link_libraries(terranni  nni ${Boost_LIBRARIES} ${EXTRA_LIBS} ${TPIE_LIBRARY})
target_link_libraries(terranni-2d  nni tpie ${Boost_LIBRARIES} ${EXTRA_LIBS} ${LIBLAS_LIBRARY} ${TPIE_LIBRARY})
target_link_libraries(terranni-2d  nni ${Boost_LIBRARIES} ${EXTRA_LIBS} ${TPIE_LIBRARY})
