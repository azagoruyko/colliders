# By Alexander Zagoruyko. No rights reserved!

if(NOT DEFINED MAYA_DEVKIT_DIR)
    set(MAYA_DEVKIT_DIR "" CACHE PATH "Maya devkit directory (where include/lib reside)")
endif()

# OS Specific environment setup
set(MAYA_COMPILE_DEFINITIONS "REQUIRE_IOSTREAM;_BOOL;_TBB_NO_IMPLICIT_LINKAGE;_TBB_MAKE_EXCEPTION_PTR_PRESENT")
set(MAYA_LIB_PREFIX "")
set(MAYA_LIB_SUFFIX ".lib")

if(WIN32)
    # Windows
    set(MAYA_COMPILE_DEFINITIONS "${MAYA_COMPILE_DEFINITIONS};NT_PLUGIN")
    set(MAYA_PLUGIN_EXTENSION ".mll")

elseif(APPLE)
    # Apple
    set(MAYA_COMPILE_DEFINITIONS "${MAYA_COMPILE_DEFINITIONS};OSMac_")
    set(MAYA_PLUGIN_EXTENSION ".bundle")
    set(MAYA_LIB_SUFFIX ".dylib")
else()
    # Linux
    set(MAYA_COMPILE_DEFINITIONS "${MAYA_COMPILE_DEFINITIONS};LINUX")
    set(MAYA_PLUGIN_EXTENSION ".so")
    set(MAYA_LIB_SUFFIX ".so")
    set(MAYA_LIB_PREFIX "lib")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Maya DEFAULT_MSG MAYA_DEVKIT_DIR)

function(MAYA_PLUGIN target)
    target_include_directories(${target} PUBLIC ${MAYA_DEVKIT_DIR}/include)

    set(MAYA_LIBS Foundation OpenMaya OpenMayaAnim OpenMayaFX OpenMayaRender OpenMayaUI clew tbb)
    foreach(MAYA_LIB ${MAYA_LIBS})
        target_link_libraries(${target} PUBLIC ${MAYA_DEVKIT_DIR}/lib/${MAYA_LIB_PREFIX}${MAYA_LIB}${MAYA_LIB_SUFFIX})
    endforeach()

    set_target_properties(${target} PROPERTIES
        COMPILE_DEFINITIONS "${MAYA_COMPILE_DEFINITIONS}"
        LINK_FLAGS "/export:initializePlugin /export:uninitializePlugin"
        PREFIX ""
        SUFFIX ${MAYA_PLUGIN_EXTENSION})
endfunction()