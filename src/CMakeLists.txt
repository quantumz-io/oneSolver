find_package(Boost 1.71 COMPONENTS program_options REQUIRED)

set(MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/src")
set(SOURCE_DIRECTORY "${MODULE_DIRECTORY}")

file(GLOB SOURCES  "${SOURCE_DIRECTORY}/*.cpp")
foreach(sourcefile ${SOURCES})
    get_filename_component(filename ${sourcefile} NAME_WE) 
    add_executable(${filename} ${sourcefile})
    target_include_directories(${filename} PRIVATE "${PROJECT_SOURCE_DIR}/include")
    target_link_libraries(${filename} ${LIBRARY_NAME})
    target_link_libraries(${filename} Boost::program_options)
endforeach(sourcefile ${SOURCES})
