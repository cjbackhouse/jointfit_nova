############   top level config
cmake_minimum_required(VERSION 3.10)
project(CAFAna)
#set(VERSION "v00.00.00")

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG")  # needed to ensure assert() causes an abort() even in optimized builds
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-undefined")

# ROOT is a non-negotiable dependency
find_package(ROOT REQUIRED)

# this is quite lazy -- it links ALL the ROOT libraries into the designated one.
# in principle we should figure out which ones we actually need and link only those.
# but.... some other time
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

# by default make the bin/ and lib/ dirs in the top of the source dir.
# can by overridden on command line if desired
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

###########   now go into the subdirs and do the actual work
add_subdirectory(OscLib)
add_subdirectory(CAFAna)