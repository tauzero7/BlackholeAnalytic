
cmake_minimum_required(VERSION 3.5.1)
project(BlackholeAnalytic)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# ------------------------
#  GNU Scientific Library
# ------------------------
set(GSL_ROOT_DIR  /usr  CACHE  FILEPATH  "GNU Scientific Library")

execute_process(COMMAND ${GSL_ROOT_DIR}/bin/gsl-config --cflags
    OUTPUT_VARIABLE GSL_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
    
execute_process(COMMAND ${GSL_ROOT_DIR}/bin/gsl-config --libs
    OUTPUT_VARIABLE GSL_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
#message(${GSL_CFLAGS} ${GSL_LIBS})

set(GSL_LIB_DIR  ${GSL_ROOT_DIR}/lib  CACHE FILEPATH "GSL library path")


# ------------------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${GSL_CFLAGS}")
include_directories(src)


# ------------------------
file(GLOB bhsrc_files  src/bhAnalytic.*)


# ------------------------
#  Ring exe
# ------------------------
add_executable(Ring main/ring.cpp ${bhsrc_files})
target_link_libraries(Ring ${GSL_LIBS})
set_target_properties(Ring  PROPERTIES LINK_FLAGS "-Wl,-rpath ${GSL_LIB_DIR}")

# ------------------------
#  Geodesics exe
# ------------------------
add_executable(Geodesics main/geodesics.cpp ${bhsrc_files})
target_link_libraries(Geodesics ${GSL_LIBS})
set_target_properties(Geodesics  PROPERTIES LINK_FLAGS "-Wl,-rpath ${GSL_LIB_DIR}")
