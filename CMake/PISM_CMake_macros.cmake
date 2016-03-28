# This file contains CMake macros used in the root CMakeLists.txt

# Set CMake variables to enable rpath
macro(pism_use_rpath)
  ## Use full RPATH, with this setting Pism libraries cannot be moved after installation
  ## but the correct libraries will always be found regardless of LD_LIBRARY_PATH
  ## in use, i.e. don't skip the full RPATH for the build tree
  set (CMAKE_SKIP_BUILD_RPATH FALSE)
  # when building, don't use the install RPATH already
  # (but later on when installing)
  set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
  # the RPATH to be used when installing
  set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
  # add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

  # Mac OS X install_name fix:
  set(CMAKE_MACOSX_RPATH 1)
  set (CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
endmacro(pism_use_rpath)

# Set CMake variables to disable rpath
macro(pism_dont_use_rpath)
  set (CMAKE_SKIP_BUILD_RPATH TRUE)
  set (CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
  set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${Pism_LIB_DIR}")
  set (CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
endmacro(pism_dont_use_rpath)

# Set CMake variables to ensure that everything is static
macro(pism_strictly_static)
  set (CMAKE_SKIP_RPATH ON CACHE BOOL "Disable RPATH completely")
  set (CMAKE_FIND_LIBRARY_SUFFIXES .a)

  set (BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared Pism libraries" FORCE)

  SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "") # get rid of -rdynamic
  SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "") # ditto

  set_property(GLOBAL PROPERTY LINK_SEARCH_END_STATIC 1)
  set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
  set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)

  pism_dont_use_rpath()
endmacro(pism_strictly_static)

# Set the revision tag if PISM was checked out using Git.
macro(pism_set_revision_tag_git)
  if (NOT Pism_VERSION)
    if (EXISTS ${Pism_SOURCE_DIR}/.git)
      find_program (GIT_EXECUTABLE git DOC "Git executable")
      mark_as_advanced(GIT_EXECUTABLE)
      execute_process (COMMAND ${GIT_EXECUTABLE} describe --always --match v?.?*
        WORKING_DIRECTORY ${Pism_SOURCE_DIR}
        OUTPUT_VARIABLE Pism_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      execute_process (COMMAND ${GIT_EXECUTABLE} --no-pager log -1 "--pretty=format:committed by %an on %ci"
        WORKING_DIRECTORY ${Pism_SOURCE_DIR}
        OUTPUT_VARIABLE Pism_COMMIT_INFO
        OUTPUT_STRIP_TRAILING_WHITESPACE)
      set(Pism_VERSION "${Pism_VERSION} ${Pism_COMMIT_INFO}")
    endif (EXISTS ${Pism_SOURCE_DIR}/.git)
  endif(NOT Pism_VERSION)
endmacro(pism_set_revision_tag_git)

# Set the PISM revision tag
macro(pism_set_revision_tag)
  # Git
  pism_set_revision_tag_git()

  # Otherwise...
  if (NOT Pism_VERSION)
    set (Pism_VERSION "no-version-control")
  endif (NOT Pism_VERSION)

  set (Pism_REVISION_TAG "${Pism_BRANCH} ${Pism_VERSION}")

  message(STATUS "Configuring PISM version '${Pism_REVISION_TAG}'")
endmacro(pism_set_revision_tag)

macro(pism_set_install_prefix)
  # Allow setting a custom install prefix using the PISM_PREFIX environment variable.
  string (LENGTH "$ENV{PISM_INSTALL_PREFIX}" INSTALL_PREFIX_LENGTH)
  if (INSTALL_PREFIX_LENGTH)
    set (CMAKE_INSTALL_PREFIX $ENV{PISM_INSTALL_PREFIX} CACHE PATH "PISM install prefix" FORCE)
    message (STATUS "Setting PISM install prefix to ${CMAKE_INSTALL_PREFIX}.")
  endif()

  # Define the directory structure.
  set (Pism_BIN_DIR "bin")
  set (Pism_LIB_DIR "lib")
  set (Pism_SHARE_DIR "share/pism")
  set (Pism_DOC_DIR "share/doc/pism")
endmacro()

# Set pedantic compiler flags
macro(pism_set_pedantic_flags)
  set (DEFAULT_PEDANTIC_FLAGS "-pedantic -Wall -Wextra -Wno-cast-qual -Wundef -Wshadow -Wpointer-arith -Wno-cast-align -Wwrite-strings -Wno-conversion -Wsign-compare -Wno-redundant-decls -Wno-inline -Wno-long-long -Wmissing-format-attribute -Wmissing-noreturn -Wpacked -Wdisabled-optimization -Wmultichar -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wendif-labels -Winvalid-pch -Wmissing-field-initializers -Wvariadic-macros -Wstrict-aliasing -funit-at-a-time")
  set (DEFAULT_PEDANTIC_CFLAGS "${DEFAULT_PEDANTIC_FLAGS} -std=c99")
  set (DEFAULT_PEDANTIC_CXXFLAGS "${DEFAULT_PEDANTIC_FLAGS} -Woverloaded-virtual")
  set (PEDANTIC_CFLAGS ${DEFAULT_PEDANTIC_CFLAGS} CACHE STRING "Compiler flags to enable pedantic warnings")
  set (PEDANTIC_CXXFLAGS ${DEFAULT_PEDANTIC_CXXFLAGS} CACHE STRING "Compiler flags to enable pedantic warnings for C++")
  mark_as_advanced (PEDANTIC_CFLAGS PEDANTIC_CXXFLAGS)
  set (CMAKE_C_FLAGS_DEBUG "-g ${PEDANTIC_CFLAGS}")
  set (CMAKE_CXX_FLAGS_DEBUG "-g ${PEDANTIC_CXXFLAGS}")
endmacro(pism_set_pedantic_flags)

# Make sure that we don't create .petscrc in $HOME, because this would affect
# all PISM runs by the current user.
macro(pism_check_build_dir_location)
  if (DEFINED ENV{HOME})
    # Don't assume that HOME env var is set.
    file (TO_CMAKE_PATH $ENV{HOME} home_dir)
    file (TO_CMAKE_PATH ${PROJECT_BINARY_DIR} build_dir)

    if (${home_dir} STREQUAL ${build_dir})
      message (FATAL_ERROR
        "\n"
        "The build directory is the same as your $HOME!\n"
        "Buiding PISM here would result in a big mess. "
        "Please create a special build directory and run cmake from there.\n")
    endif()
  endif()
endmacro()

macro(pism_find_prerequisites)
  # PETSc
  find_package (PETSc)
  if (DEFINED PETSC_VERSION)
    # FindPETSc.cmake does not put PETSC_VERSION into the CMake cache,
    # so we save it here.
    set(Pism_PETSC_VERSION ${PETSC_VERSION} CACHE STRING "PETSc version")
    mark_as_advanced(Pism_PETSC_VERSION)
  endif()

  if ((DEFINED PETSC_VERSION) AND (PETSC_VERSION VERSION_LESS 3.3))
    # Force PISM to look for PETSc again if the version we just found
    # is too old:
    set(PETSC_CURRENT "OFF" CACHE BOOL "" FORCE)
    # Stop with an error message.
    message(FATAL_ERROR "\nPISM requires PETSc version 3.3 or newer (found ${PETSC_VERSION}).\n\n")
  endif()

  if ((DEFINED PETSC_VERSION) AND (PETSC_VERSION VERSION_EQUAL 3.6.0))
    # Force PISM to look for PETSc again if the version we just found
    # is not supported
    set(PETSC_CURRENT "OFF" CACHE BOOL "" FORCE)
    # Stop with an error message.
    message(FATAL_ERROR "\nPISM does not support PETSc ${PETSC_VERSION}. Please install PETSc <= 3.5.4 or PETSc > 3.6.0.\n\n")
  endif()

  # MPI
  # Use the PETSc compiler as a hint when looking for an MPI compiler
  # FindMPI.cmake changed between 2.8.4 and 2.8.5, so we try to support both...
  if (${CMAKE_VERSION} VERSION_LESS "2.8.5")
    set (MPI_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")
    find_package (MPI REQUIRED)

    set (MPI_C_INCLUDE_PATH "${MPI_INCLUDE_PATH}" CACHE STRING "MPI include directories (semicolon-separated list)")
    set (MPI_C_LIBRARIES "${MPI_LIBRARY};${MPI_EXTRA_LIBRARY}" CACHE STRING "MPI libraries (semicolon-separated list)")
    mark_as_advanced(MPI_C_INCLUDE_PATH MPI_C_LIBRARIES)
    message (STATUS
      "Note: Please upgrade CMake to version 2.8.5 or later if the build fails with undefined symbols related to MPI.")
  else ()
    set (MPI_C_COMPILER ${PETSC_COMPILER} CACHE FILEPATH "MPI compiler. Used only to detect MPI compilation flags.")
    find_package (MPI REQUIRED)
  endif()

  # Other required libraries
  find_package (UDUNITS2 REQUIRED)
  find_package (GSL REQUIRED)
  find_package (NetCDF REQUIRED)

  # Optional libraries
  if (Pism_USE_PNETCDF)
    find_package (PNetCDF)
  endif()
  if (Pism_USE_PARALLEL_HDF5)
    find_package (HDF5 COMPONENTS C HL)
  endif()
  find_package (FFTW REQUIRED)    # NOT optional?
  if (Pism_USE_PROJ4)
    find_package (PROJ4)
  endif()

  # Use TAO included in PETSc 3.5.
  if (Pism_PETSC_VERSION VERSION_LESS "3.5")
    message(STATUS "Disabling TAO-based inversion tools. Install PETSc 3.5 or later to use them.")
  else()
    set (Pism_USE_TAO ON CACHE BOOL "Use TAO in inverse solvers.")
    if (Pism_USE_TAO)
      message(STATUS "PETSc 3.5 and later includes TAO; using it...")
    else()
      message(STATUS "Pism_USE_TAO is OFF. Inverse solvers using the TAO library will not be built.")
    endif()
  endif()

  # Try to find netcdf_par.h. We assume that NetCDF was compiled with
  # parallel I/O if this header is present.
  find_file(NETCDF_PAR_H netcdf_par.h HINTS ${NETCDF_INCLUDES} NO_DEFAULT_PATH)

  # Set default values for build options
  if (NOT NETCDF_PAR_H)
    set (Pism_USE_PARALLEL_NETCDF4 OFF CACHE BOOL "Enables parallel NetCDF-4 I/O." FORCE)
    message(STATUS "Selected NetCDF library does not support parallel I/O.")
  endif()

  if (NOT PNETCDF_FOUND)
    set (Pism_USE_PNETCDF OFF CACHE BOOL "Enables parallel NetCDF-3 I/O using PnetCDF." FORCE)
  endif()

  if (NOT HDF5_FOUND)
    set (Pism_USE_PARALLEL_HDF5 OFF CACHE BOOL "Enables parallel HDF5 I/O." FORCE)
  elseif(NOT HDF5_IS_PARALLEL)
    set (Pism_USE_PARALLEL_HDF5 OFF CACHE BOOL "Enables parallel HDF5 I/O." FORCE)
    message (STATUS "Selected HDF5 library does not support parallel I/O.")
  endif()

  if (PROJ4_FOUND)
    set (Pism_USE_PROJ4 ON CACHE BOOL
      "Use Proj.4 to compute cell areas, longitude, and latitude.")
  endif()

endmacro()

macro(pism_set_dependencies)

  # Set include and library directories for *required* libraries.
  include_directories (
    ${PETSC_INCLUDES}
    ${FFTW_INCLUDE_DIRS}
    ${FFTW_INCLUDES}
    ${GSL_INCLUDES}
    ${UDUNITS2_INCLUDES}
    ${NETCDF_INCLUDES}
    ${MPI_C_INCLUDE_PATH})

  # Use option values to set compiler and linker flags
  set (Pism_EXTERNAL_LIBS "")

  # required libraries
  list (APPEND Pism_EXTERNAL_LIBS
    ${PETSC_LIBRARIES}
    ${UDUNITS2_LIBRARIES}
    ${FFTW_LIBRARIES}
    ${GSL_LIBRARIES}
    ${NETCDF_LIBRARIES}
    ${MPI_C_LIBRARIES})

  # optional libraries
  if (Pism_USE_PROJ4)
    include_directories (${PROJ4_INCLUDES})
    list (APPEND Pism_EXTERNAL_LIBS ${PROJ4_LIBRARIES})
  endif()

  if (Pism_USE_PNETCDF)
    include_directories (${PNETCDF_INCLUDES})
    list (APPEND Pism_EXTERNAL_LIBS ${PNETCDF_LIBRARIES})
  endif()

  # Put HDF5 includes near the beginning of the list. (It is possible that the system has
  # more than one HDF5 library installed--- one serial, built with NetCDF, and one parallel.
  # We want to use the latter.)
  if (Pism_USE_PARALLEL_HDF5)
    include_directories (BEFORE ${HDF5_C_INCLUDE_DIR})
    list (APPEND Pism_EXTERNAL_LIBS ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  endif()

  # Hide distracting CMake variables
  mark_as_advanced(file_cmd MPI_LIBRARY MPI_EXTRA_LIBRARY
    CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT
    MAKE_EXECUTABLE HDF5_DIR NETCDF_PAR_H)

endmacro()

include(CheckCXXSourceCompiles)

# Check if shared_ptr is in <memory> as std::shared_ptr. If it is not,
# assume that we have to use <tr1/memory> and std::tr1::shared_ptr.
macro(pism_check_shared_ptr)
  set(SHARED_PTR_TEST_SRC "
#include <memory>

int main(int argc, char **argv) {
  std::shared_ptr<double> shared_double(new double);
  return 0;
}
")
  check_cxx_source_compiles("${SHARED_PTR_TEST_SRC}" PISM_SHARED_PTR)
  if (PISM_SHARED_PTR)
    set(Pism_USE_TR1 OFF CACHE BOOL "If 'ON', #include <tr1/memory>, otherwise #include <memory>.")
  else()
    set(Pism_USE_TR1 ON CACHE BOOL "If 'ON', #include <tr1/memory>, otherwise #include <memory>." FORCE)
  endif()
endmacro()
