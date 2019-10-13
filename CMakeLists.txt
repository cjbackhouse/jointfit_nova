############   top level config
cmake_minimum_required(VERSION 3.10)
project(CAFAna)
#set(VERSION "v00.00.00")

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -NDEBUG")  # needed to ensure assert() causes an abort() even in optimized builds

# ROOT is a non-negotiable dependency
find_package(ROOT REQUIRED)

function(link_root target)
    unset(libs CACHE)
    list(APPEND libs ${ROOT_LIBRARIES})
    foreach(lib ${libs})
        target_link_libraries(${target} PUBLIC ${lib})
    endforeach()
endfunction()


set(TARGET_LIBDIR lib)
set(TARGET_BINDIR bin)

###########   directory configuration

include_directories(
        ${CMAKE_SOURCE_DIR}
        ${ROOT_INCLUDE_DIRS}
)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

###########   now go into the subdirs and do the actual work
add_subdirectory(OscLib)
add_subdirectory(CAFAna)