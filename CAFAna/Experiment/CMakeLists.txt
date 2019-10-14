set(LIBRARY CAFAnaExperiment)

file(GLOB HEADER_FILES *.h)
file(GLOB SOURCES *.cxx)

add_library(${LIBRARY} SHARED
        ${HEADER_FILES}
        ${SOURCES}
)
link_root(${LIBRARY})

target_link_libraries(${LIBRARY} PUBLIC CAFAnaCore)

install(TARGETS ${LIBRARY} LIBRARY DESTINATION ${TARGET_LIBDIR})
install_headers("${HEADER_FILES}")
