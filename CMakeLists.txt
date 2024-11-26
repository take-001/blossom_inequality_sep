cmake_minimum_required(VERSION 3.10)
project(BlossomSeparation)

set(CMAKE_C_STANDARD 11)

# Source files for the library
set(SOURCES
    blossom_detector.c
    src/blossom.c
    src/cliqwork.c
    src/skeleton.c
    src/cutpool.c
    src/allocrus.c
    src/urandom.c
    src/gomoryhu.c
    src/sortrus.c
    src/safe_io.c
    src/cliqhash.c
    src/genhash.c
    src/cut_st.c
    src/util.c
)

# Create the shared library
add_library(blossom_separation SHARED ${SOURCES})

# Include headers
target_include_directories(blossom_separation PRIVATE ${CMAKE_SOURCE_DIR}/INCLUDE)

# Link with math library
target_link_libraries(blossom_separation PRIVATE m)

# Add compiler options for shared library
set_target_properties(blossom_separation PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Set the output directory for the shared library
set_target_properties(blossom_separation PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)


add_executable(blossom_detector ${SOURCES})

target_include_directories(blossom_detector PRIVATE ${CMAKE_SOURCE_DIR}/INCLUDE)

target_link_libraries(blossom_detector PRIVATE m)
