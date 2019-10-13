set(LIBRARY OscLib)

file(GLOB HEADER_FILES *.h)
file(GLOB SOURCES *.cxx)

add_library(${LIBRARY} SHARED
        ${HEADER_FILES}
        ${SOURCES}
        )

install(TARGETS ${LIBRARY} LIBRARY DESTINATION ${TARGET_LIBDIR})
