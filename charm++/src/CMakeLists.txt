set(CMAKE_CXX_COMPILER "${CHARMC}")
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB TABASCO_SOURCES *.cc)
add_executable(tabasco ${TABASCO_SOURCES})
add_dependencies(tabasco charmc)

file(GLOB TABASCO_CHARM_MODULES *.ci)
foreach(_CI ${TABASCO_CHARM_MODULES})
	get_filename_component(_TRUNC "${_CI}" NAME_WE)
	add_custom_command(OUTPUT ${_TRUNC}.decl.h ${_TRUNC}.def.h 
		COMMAND ${CHARMC} ${_CI}
		COMMENT "Compiling ${_CI}"
		DEPENDS charmc)
	add_custom_target(${_TRUNC} DEPENDS ${_TRUNC}.decl.h
		${_TRUNC}.def.h)
	add_dependencies(tabasco ${_TRUNC})
endforeach()
