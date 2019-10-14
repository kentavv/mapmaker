cmake_minimum_required(VERSION 3.11)
project(mapmaker C)

#set(CMAKE_C_STANDARD 17)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

SET(CMAKE_C_FLAGS "-g -Wall -Wextra -Wno-unused-parameter -Wno-clobbered -Wno-format-overflow -pedantic -std=gnu17")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

set(CMAKE_VERBOSE_MAKEFILE ON)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#OPTION(USE_OpenMP "Use OpenMP" ON)
#IF (USE_OpenMP)
#    FIND_PACKAGE(OpenMP)
#    IF (OPENMP_FOUND)
#        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    ENDIF ()
#ENDIF ()

include_directories(. lib)
add_compile_options(-Wall -Wextra -Wno-sign-compare)

add_library(mylib STATIC
        lib/eqn.c
        lib/eqn.h
        lib/iolib.c
        lib/iolib.h
        lib/mathlib.c
        lib/mathlib.h
        lib/memlib.c
        lib/memlib.h
        lib/misclib.h
        lib/msglib.c
        lib/msglib.h
        lib/shell.c
        lib/shell.h
        lib/stats.c
        lib/stats.h
        lib/strlib.c
        lib/strlib.h
        lib/syscode.c
        lib/syscode.h
        lib/system.h
        lib/table.c
        lib/table.h
        )

add_executable(
        makehelp
        makehelp/makehelp.c
)

add_executable(
        mapmaker
        mapm/auto_cmd.c
        mapm/chroms.c
        mapm/ctm.c
        mapm/database.c
        mapm/info.c
        mapm/lowlevel.h
        mapm/main.c
        mapm/map_info.h
        mapm/mapm.h
        mapm/mapmhelp.c
        mapm/maps.c
        mapm/npt_cmds.c
        mapm/orders.c
        mapm/print.c
        mapm/ps_maps.c
        mapm/quick23.c
        mapm/reader.c
        mapm/sequence.c
        mapm/state.c
        mapm/sys_cmds.c
        mapm/toplevel.h
        mapm/two_cmds.c
)

add_executable(
        qtl
        quant/qcmds.c
        quant/qcontext.c
        quant/qctm.c
        quant/qdata.c
        quant/qlow.h
        quant/qmap.h
        quant/qprint.c
        quant/qps_scan.c
        quant/qraw.c
        quant/qseq.c
        quant/qtl.h
        quant/qtlhelp.c
        quant/qtop.c
        quant/qtop.h
        quant/qwiggle.c
)

set(mapmaker_help_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/mapm/mapmhelp.src)
set(qtl_help_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/quant/qtlhelp.src)

add_custom_target(build_mapmaker_help
        COMMAND makehelp ./mapm/mapmhelp.src ./mapm/mapmhelp.c ./mapmaker.help ./mapmaker.ps
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        SOURCES ${mapmaker_help_SOURCES}
        )

add_custom_target(build_qtl_help
        COMMAND makehelp ./quant/qtlhelp.src ./quant/qtlhelp.c ./qtl.help ./qtl.ps
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        SOURCES ${mapmaker_help_SOURCES}
        )

add_dependencies(mapmaker build_mapmaker_help)
add_dependencies(qtl build_qtl_help)

target_link_libraries(makehelp mylib m readline termcap)
target_link_libraries(mapmaker mylib m readline termcap)
target_link_libraries(qtl mylib m readline termcap)
