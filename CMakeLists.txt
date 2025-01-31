# Specify the minimum version of CMake
cmake_minimum_required(VERSION 3.15)

# Project name and language
project(Costo LANGUAGES Fortran C)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set (CMAKE_Fortran_FLAGS "-Wall -O0")
elseif (Fortran_COMPILER_NAME MATCHES "ifx.*")
  set (CMAKE_Fortran_FLAGS "-heap-arrays -g -O0")
endif()


# Specify the Fortran standard
#set(CMAKE_Fortran_STANDARD 90)
#set(CMAKE_Fortran_STANDARD_REQUIRED TRUE)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED TRUE)



# Define the source directories
set(MODULE_I_DIR "${CMAKE_SOURCE_DIR}/src/INTERF")
set(MODULE_DECLA_DIR "${CMAKE_SOURCE_DIR}/src/DECLA")
set(MODULE_NRF_DIR "${CMAKE_SOURCE_DIR}/src/NRF90")
set(PROCEDURE_DIR "${CMAKE_SOURCE_DIR}/src/procedures")
set(C_HELPER_DIR "${CMAKE_SOURCE_DIR}/src/lib_C")
set(MAIN_FILE "${CMAKE_SOURCE_DIR}/src/joint_inv.f90")
set(BIN "${CMAKE_SOURCE_DIR}/bin")

set(NRF_S "${MODULE_NRF_DIR}/mynrutil.f90" "${MODULE_NRF_DIR}/nrtype.f90" "${MODULE_NRF_DIR}/nrutil.f90")

# Collect all source files
file(GLOB MODULE_I_SOURCES "${MODULE_I_DIR}/*.f90")
file(GLOB MODULE_DECLA_SOURCES "${MODULE_DECLA_DIR}/*.f")
#file(GLOB MODULE_NRF_SOURCES "${MODULE_NRF_DIR}/*.f90")
file(GLOB PROCEDURE_SOURCES "${PROCEDURE_DIR}/*.f90")
file(GLOB C_SOURCES "${C_HELPER_DIR}/*.c")


set(SOURCES ${MAIN_FILE} ${MODULE_DECLA_SOURCES}  ${MODULE_I_SOURCES}  ${NRF_S}   ${PROCEDURE_SOURCES}  ${C_SOURCES})

#add_executable(costo ${MAIN_FILE})
add_executable(costo ${SOURCES})

# Specify module output directory
set_target_properties(costo PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")

# Compiler options
target_compile_options(costo PRIVATE -O0)

# Optional: Debug or Release configurations
if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_options(costo PRIVATE -g)
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(costo PRIVATE -O3)
endif()
