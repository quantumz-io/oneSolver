find_package(Boost 1.71 COMPONENTS program_options REQUIRED)

set(APP_MODULE_DIRECTORY "${PROJECT_SOURCE_DIR}/app")
set(APP_SOURCE_DIRECTORY "${APP_MODULE_DIRECTORY}")

file(GLOB APP_SOURCES  "${APP_SOURCE_DIRECTORY}/*.cpp")
foreach(sourcefile ${APP_SOURCES})
    get_filename_component(filename ${sourcefile} NAME_WE) 
    add_executable(${filename} ${sourcefile})
    target_include_directories(${filename} PRIVATE "${PROJECT_SOURCE_DIR}/include")
    target_link_libraries(${filename} ${LIBRARY_NAME})
    target_link_libraries(${filename} Boost::program_options)
endforeach(sourcefile ${APP_SOURCES})
