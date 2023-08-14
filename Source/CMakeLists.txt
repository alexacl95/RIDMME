set(RIDMME_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/Private/AutomaticDifferentiation.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Private/MatrixRIDMME.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Private/Metrics.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Private/Montecarlo.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Private/NumericalMethod.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Private/Refinement.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/Private/RIDMME_Interval_Autodiff.cpp")
set(RIDMME_HEADERS
    "${CMAKE_CURRENT_SOURCE_DIR}/Public/AutomaticDifferentiation.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/Public/MatrixRIDMME.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/Public/Metrics.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/Public/Montecarlo.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/Public/NumericalMethod.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/Public/Refinement.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/Public/RIDMME_Interval_Autodiff.h")

set(RIDMME_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/Public/")

# MyLib Library

add_library(${RIDMME_LIBRARY_NAME} STATIC ${RIDMME_SOURCES} ${RIDMME_HEADERS})
target_include_directories(
    ${RIDMME_LIBRARY_NAME}
    PUBLIC ${RIDMME_INCLUDES})

if(OpenMP_CXX_FOUND)
target_link_libraries(${RIDMME_LIBRARY_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(
    ${RIDMME_LIBRARY_NAME}
    PUBLIC bias::bias)
target_set_warnings(
    TARGET
    ${RIDMME_LIBRARY_NAME}
    ENABLE
    ${ENABLE_WARNINGS}
    AS_ERRORS
    ${ENABLE_WARNINGS_AS_ERRORS})

