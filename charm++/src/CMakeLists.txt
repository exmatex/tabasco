file(GLOB TABASCO_CHARM_MODULES *.ci)
foreach(_CI ${TABASCO_CHARM_MODULES})
	get_filename_component(_TRUNC "${_CI}" NAME_WE)
	message("${_CI} ${_TRUNC}")
	add_custom_command(OUTPUT ${_TRUNC}.decl.h ${_TRUNC}.def.h 
		COMMAND ${CHARMC} ${_CI}
		COMMENT "Compiling ${_CI}")
	add_custom_target(${_TRUNC} ALL DEPENDS ${_TRUNC}.decl.h
		${_TRUNC}.def.h)
endforeach()
set(CMAKE_CXX_COMPILER "${CHARMC}")
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB TABASCO_SOURCES *.cc)
add_executable(tabasco ${TABASCO_SOURCES})
add_dependencies(tabasco charmc)
