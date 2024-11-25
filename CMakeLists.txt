cmake_minimum_required(VERSION 3.10)
project(BlossomDetector)


set(CMAKE_C_STANDARD 11)


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



add_executable(blossom_detector ${SOURCES})

target_include_directories(blossom_detector PRIVATE ${CMAKE_SOURCE_DIR}/INCLUDE)

target_link_libraries(blossom_detector PRIVATE m)

