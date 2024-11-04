cmake_minimum_required(VERSION 3.16.3)

project(libnabo)

if("$ENV{ROS_VERSION}" STREQUAL "1")

# Source and header files.
set(NABO_SRC
  nabo/nabo.cpp
  nabo/brute_force_cpu.cpp
  nabo/kdtree_cpu.cpp
  nabo/kdtree_opencl.cpp
)

file(GLOB_RECURSE NABO_HEADERS "nabo/*.h")

# Build type options:
#	* catkin.
#	* cmake.
#	* ament (future).
set(BUILD_TYPE "catkin")
if (BUILD_TYPE STREQUAL "catkin")
  include(CatkinBuild.cmake)
else()
  # For a pure-cmake build, follow the instructions on the package.xml file, which
  # specify how to set the buildtool to cmake and export the package accordingly.
  include(PureCMakeBuild.cmake)
endif()

else() # ROS version 2

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall -Wextra -Wpedantic -Werror=return-type)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(PACKAGE_DEPENDENCIES
  Eigen3
  )

find_package(ament_cmake REQUIRED)

foreach(PKG ${PACKAGE_DEPENDENCIES})
  find_package(${PKG} REQUIRED)
endforeach()

if(USE_OPENMP)
  message(FATAL_ERROR "We are not going to build with openmp on colcon.")
endif(USE_OPENMP)

########################
## Library definition ##
########################

add_library(nabo
  nabo/nabo.cpp
  nabo/brute_force_cpu.cpp
  nabo/kdtree_cpu.cpp
  )

target_include_directories(nabo
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    $<INSTALL_INTERFACE:include/>
  )

target_link_libraries(nabo
  Eigen3::Eigen
  )

# NOTE(apoghosov): it seems that at the time of porting to ROS2 we have the following setup in our build pipeline:
# Release - -O3 by default
# RelWithDebInfo - -O2 by default
# but I want to follow what is written in the catkin build variant ... -O3 everywhere except debug. Note that
# specifiying it here will overide the previous option ... last specified is applied
target_compile_options(nabo PRIVATE
  "$<$<NOT:$<CONFIG:Debug>>:-O3>"
  )

set_target_properties(nabo
  PROPERTIES
    POSITION_INDEPENDENT_CODE ON
  )

ament_target_dependencies(nabo ${PACKAGE_DEPENDENCIES})

#############
## Install ##
#############

install(DIRECTORY nabo/
  DESTINATION include/nabo
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "opencl" EXCLUDE
  )

install(TARGETS nabo
  EXPORT export_${PROJECT_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  )

ament_export_targets(export_${PROJECT_NAME} HAS_LIBRARY_TARGET)
ament_export_dependencies(${PACKAGE_DEPENDENCIES})

#############
## Testing ##
#############
  
if(BUILD_TESTING)
  
  find_package(ament_cmake_gtest REQUIRED)
  
  ament_add_gtest(test_${PROJECT_NAME}
    tests/empty_test.cpp
    )

  target_link_libraries(test_${PROJECT_NAME} nabo)

  ###################
  ## Code_coverage ##
  ###################

  find_package(cmake_code_coverage QUIET)

  if(cmake_code_coverage_FOUND)
    add_gtest_coverage(TEST_BUILD_TARGETS test_${PROJECT_NAME})
  endif()

endif()

ament_package()

endif()
