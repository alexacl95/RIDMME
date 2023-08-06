# Project/Library Names
set(EXECUTABLE_NAME "example")
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/Source/Public")

set(SAMPLE_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Models.cpp")
set(SAMPLE_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/Models.h")

# MyLib Library
add_executable(${EXECUTABLE_NAME} ${SAMPLE_SOURCES} ${SAMPLE_HEADERS})

target_link_libraries(
    ${EXECUTABLE_NAME}
    PUBLIC ${RIDMME_LIBRARY_NAME})
target_set_warnings(
    TARGET
    ${EXECUTABLE_NAME}
    ENABLE
    ${ENABLE_WARNINGS}
    AS_ERRORS
    ${ENABLE_WARNINGS_AS_ERRORS})

