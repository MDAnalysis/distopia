# distopia python layer
add_cython_target(pydistopia CXX PY3)
add_library(pydistopia MODULE ${pydistopia})
python_extension_module(pydistopia)
target_link_libraries(pydistopia libdistopia)

if(APPLE)
    set_target_properties(pydistopia PROPERTIES INSTALL_RPATH "@loader_path")
else()
    set_target_properties(pydistopia PROPERTIES INSTALL_RPATH "\$ORIGIN")
endif()

install(TARGETS pydistopia LIBRARY DESTINATION distopia)