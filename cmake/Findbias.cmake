#
# Helper function for finding BIAS.

#
set(_bias_version "22.0.8")
message("-- Looking for BIAS version: ${_bias_version}")

set(_bias_approot "${PROJECT_SOURCE_DIR}/third-party/Profil-2.0.8")
set(_bias_libdir "lib")
set(_bias_libname "libProfilPackages.a")
set(BIAS_DEPENDENCIES "libProfil.a" "libBias.a" "liblr.a")


set(BIAS_DEPENDENCIES_LIBS)

set(bias_LINKER_FLAGS_LIST
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,SHARED_LIBRARY>:>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,MODULE_LIBRARY>:>"
    "$<$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>:>")

# should point the the BIAS SDK installation dir
set(_bias_root "${_bias_approot}")
# find header dir and libs
find_path(
    bias_INCLUDE_DIR "Interval.h"
    PATHS ${_bias_root}
    PATH_SUFFIXES "include")

find_library(
    bias_LIBRARY ${_bias_libname}
    PATHS ${_bias_root}
    PATH_SUFFIXES ${_bias_libdir})
    
foreach(_DEP ${BIAS_DEPENDENCIES})
    find_library(
        FOUND_DEP ${_DEP}
        PATHS ${_bias_root}
        PATH_SUFFIXES ${_bias_libdir}
        NO_CACHE)
    set(BIAS_DEPENDENCIES_LIBS ${BIAS_DEPENDENCIES_LIBS} ${FOUND_DEP})
    unset(FOUND_DEP)
endforeach()

#if(NOT BIAS_DEPENDENCIES_LIBS)
 #   message(SEND_ERROR "not found dependencies")
#endif()

if(bias_INCLUDE_DIR AND bias_LIBRARY)
    set(bias_FOUND YES)
    message("-- found bias: ${_bias_version} (found version \"${_bias_version}\")")

    set(bias_LINK_LIBS ${bias_LIBRARY} ${BIAS_DEPENDENCIES_LIBS})

    if(NOT
       ${CMAKE_VERSION}
       VERSION_LESS
       "3.0")
        # Target approach
        if(NOT TARGET bias::bias)
            add_library(bias::bias INTERFACE IMPORTED)
            if(bias_INCLUDE_DIR)
                set_target_properties(
                    bias::bias PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                              "${bias_INCLUDE_DIR}")
            endif()
            set_property(
                TARGET bias::bias
                PROPERTY INTERFACE_LINK_LIBRARIES
                         "${bias_LINK_LIBS};${bias_LINKER_FLAGS_LIST}")

        endif()
    endif()
else()
    set(bias_FOUND NO)
    message(SEND_ERROR "bias not found")
endif()
