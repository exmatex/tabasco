set(CMAKE_CXX_COMPILER "${CHARMC}")
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB TABASCO_SOURCES *.cc)
include_directories($ENV{BOOST_INCLUDES})
add_executable(tabasco ${TABASCO_SOURCES})

file(GLOB TABASCO_CHARM_MODULES *.ci)
foreach(_CI ${TABASCO_CHARM_MODULES})
	get_filename_component(_TRUNC "${_CI}" NAME_WE)
	add_custom_command(OUTPUT ${_TRUNC}.decl.h ${_TRUNC}.def.h 
		COMMAND ${CHARMC} ${_CI}
		COMMENT "Compiling ${_CI}")
	add_custom_target(${_TRUNC} DEPENDS ${_TRUNC}.decl.h
		${_TRUNC}.def.h)
	add_dependencies(tabasco ${_TRUNC})
endforeach()
include_directories(${CMAKE_SOURCE_DIR}/CoEVP)
include_directories(${CMAKE_SOURCE_DIR}/CoEVP/CM/include)
target_link_libraries(tabasco lulesh_core cm)
