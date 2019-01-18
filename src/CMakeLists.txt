set (SOURCES_LIB
    FreestreamMilne.cpp    
    )

set (SOURCES
    RunWrapper.cpp
    )

add_library(freestream_lib SHARED ${SOURCES_LIB})
add_executable(RunWrapper ${SOURCES})
target_link_libraries (freestream_lib ${LIBS})
target_link_libraries (RunWrapper ${LIBS})

install(TARGETS RunWrapper DESTINATION ${CMAKE_HOME_DIRECTORY})
