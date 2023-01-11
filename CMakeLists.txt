cmake_minimum_required(VERSION 3.13)

project(LAPACK-HOMEWORK)


option(BUILD_LAPACK_FROM_SRC "Build Full LAPACK(lapack, blas, cblas and lapcke)" OFF)
option(BUILD_CBLAS_FROM_LAPACK_SRC "Build Cblas from LAPACK source code" OFF)
option(BUILD_LAPACKE_FROM_LAPACK_SRC "Build LAPACKE from LAPACK source code" OFF)

option(BUILD_GOTO2 "Build with GotoBLAS2" OFF)
option(BUILD_ATLAS "Build with ATLAS" OFF)

link_libraries(boost_program_options)
add_executable(main main.cc)
target_compile_options(main PUBLIC -std=c++17)
list(APPEND LAPACKE_TARGETS main)


if (BUILD_LAPACK_FROM_SRC)
    message(STATUS "Force to build full LAPACK from LAPACKE source")

    list(LENGTH CBLAS_TARGETS CBLAS_LIST_LEN)
    if(NOT (CBLAS_LIST_LEN EQUAL 0))
        set(CBLAS ON)
    endif()

    list(LENGTH LAPACKE_TARGETS LAPACKE_LIST_LEN)
    if(NOT (LAPACKE_LIST_LEN EQUAL 0))
        set(CBLAS ON)
        set(LAPACKE ON)
    endif()
elseif (BUILD_GOTO2)
    find_library(LAPACK_EXISTS lapack REQUIRED)
    message(STATUS "LAPACK lib found: " ${LAPACK_EXISTS})
elseif (BUILD_ATLAS)
    find_library(CBLAS_EXISTS cblas REQUIRED)
    message(STATUS "ATLAS CBLAS lib found: " ${CBLAS_EXISTS})
    find_library(ATLAS_EXISTS atlas REQUIRED)
    message(STATUS "ATLAS lib found: " ${ATLAS_EXISTS})
    find_library(ATLAS_LAPACK_EXISTS lapack_atlas REQUIRED)
    message(STATUS "ATLAS LAPACK lib found: " ${ATLAS_LAPACK_EXISTS})
else()

    # By default build from lapack
    find_library(BLAS_EXISTS blas)
    if (${BLAS_EXISTS} STREQUAL "BLAS_EXISTS-NOTFOUND")
        message(STATUS "BLAS not found")
        set(BUILD_LAPACK_FROM_SRC ON)
    endif()

    find_library(LAPACK_EXISTS lapack)
    if (${LAPACK_EXISTS} STREQUAL "LAPACK_EXISTS-NOTFOUND")
        message(STATUS "LAPACK not found")
        set(BUILD_LAPACK_FROM_SRC ON)
    endif()

    list(LENGTH CBLAS_TARGETS CBLAS_LIST_LEN)
    if(NOT (CBLAS_LIST_LEN EQUAL 0))
        if (BUILD_CBLAS_FROM_LAPACK_SRC)
            message(STATUS "Force to build CBLAS from LAPACKE source")
            set(CBLAS ON)
            set(BUILD_LAPACK_FROM_SRC ON)
        else()
            find_library(CBLAS_EXISTS cblas)
            if (${CBLAS_EXISTS} STREQUAL "CBLAS_EXISTS-NOTFOUND")
                message(STATUS "CBLAS not found")
                set(CBLAS ON)
                set(BUILD_LAPACK_FROM_SRC ON)
            else()
                message(STATUS "CBLAS lib found: " ${CBLAS_EXISTS})
            endif()
        endif()
    endif()
    

    list(LENGTH LAPACKE_TARGETS LAPACKE_LIST_LEN)
    if(NOT (LAPACKE_LIST_LEN EQUAL 0))
        if (BUILD_LAPACKE_FROM_LAPACK_SRC)
            message(STATUS "Force to build LAPACKE from LAPACKE source")
            set(LAPACKE ON)
            set(BUILD_LAPACK_FROM_SRC ON)
        else()
            find_library(LAPACKE_EXISTS lapacke)
            if (${LAPACKE_EXISTS} STREQUAL "LAPACKE_EXISTS-NOTFOUND")
                message(STATUS "LAPACKE not found")
                set(LAPACKE ON)
                set(BUILD_LAPACK_FROM_SRC ON)
            else()
                message(STATUS "LAPACKE lib found: " ${LAPACKE_EXISTS})
            endif()
        endif()
    endif()
    


    
endif()


set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

if(BUILD_LAPACK_FROM_SRC)
    message(STATUS "Enable LAPACK Source Build")
    set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
    add_subdirectory(lapack-3.11.0)

    if(CBLAS)
        foreach(target ${CBLAS_TARGETS})
            target_compile_definitions(${target} PUBLIC BUILD_LAPACK=)
            include_directories(${target} PUBLIC lapack-3.11.0/CBLAS/include)
            add_dependencies(${target} cblas)
        endforeach()
    endif()
    

    if (LAPACKE)
        foreach(target ${LAPACKE_TARGETS})
            target_compile_definitions(${target} PUBLIC BUILD_LAPACK=)
            include_directories(${target} PUBLIC lapack-3.11.0/CBLAS/include)
            include_directories(${target} PUBLIC lapack-3.11.0/LAPACKE/include)
            add_dependencies(${target} lapacke cblas)
        endforeach()
    endif()
    foreach(target ${CBLAS_TARGETS})
        target_link_libraries(${target} cblas m gfortran)
    endforeach()

    foreach(target ${LAPACKE_TARGETS})
        target_link_libraries(${target} lapacke cblas m gfortran)
    endforeach()

elseif (BUILD_GOTO2)
    message(STATUS "Build with GOTO2")
    foreach(target ${CBLAS_TARGETS})
        target_compile_definitions(${target} PRIVATE BUILD_GOTO2=)
        target_include_directories(${target} PUBLIC GotoBLAS2)
        target_link_directories(${target} PUBLIC GotoBLAS2)
        target_link_libraries(${target} goto2 m gfortran)
    endforeach()

    foreach(target ${LAPACKE_TARGETS})
        target_compile_definitions(${target} PRIVATE BUILD_GOTO2=)
        target_include_directories(${target} PUBLIC GotoBLAS2)
        target_link_directories(${target} PUBLIC GotoBLAS2)
        target_link_libraries(${target} goto2 lapack m gfortran)
    endforeach()

elseif (BUILD_ATLAS)
    message(STATUS "Build with ATLAS")
    foreach(target ${CBLAS_TARGETS})
        target_compile_definitions(${target} PUBLIC BUILD_ATLAS=)
        target_link_libraries(${target} cblas atlas lapack_atlas m gfortran)
    endforeach()    
    
    foreach(target ${LAPACKE_TARGETS})
        target_compile_definitions(${target} PUBLIC BUILD_ATLAS=)
        target_link_libraries(${target} cblas atlas lapack_atlas m gfortran)
    endforeach()
endif()





