CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(seqhax C)

# Cmake options
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
enable_testing()

# git versioning
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/libqes/cmake-modules")

INCLUDE(GitSemVer)
GetGitSemVer(VERSION)
MESSAGE(STATUS "Building seqhax version:  ${VERSION}")

# CC options
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 -DSHX_VERSION=\\\"${VERSION}\\\"")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -ggdb -Wall")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O4 -Werror")

###############################
## Find Packages and Headers ##
###############################

FIND_PACKAGE(ZLIB 1.2.5 REQUIRED)
SET(SEQHAX_DEPENDS_CFLAGS
	${SEQHAX_DEPENDS_CFLAGS} ${ZLIB_CFLAGS})
SET(SEQHAX_LIB_DEPENDS
    ${SEQHAX_LIB_DEPENDS} qes ${ZLIB_LIBRARIES})
SET(SEQHAX_DEPENDS_INCLUDE_DIRS
    ${SEQHAX_DEPENDS_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${SEQHAX_DEPENDS_INCLUDE_DIRS})

##########################
## Set Compiler Options ##
##########################

IF (CMAKE_COMPILER_IS_GNUCC)
	SET(SEQHAXWRN "${SEQHAXWRN} -Woverride-init -Wnormalized=id -Wlogical-op")
	EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion
			OUTPUT_VARIABLE GCC_VERSION)
	IF (GCC_VERSION VERSION_GREATER 4.9 OR GCC_VERSION VERSION_EQUAL 4.9)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")
	ENDIF()
ENDIF()

# Set CFLAGS
SET(SEQHAXWRN "${SEQHAXWRN} -fstack-protector-all -Wstack-protector -Wfloat-equal")
SET(SEQHAXWRN "${SEQHAXWRN} -Wundef -Wpointer-arith -Wstrict-prototypes")
SET(SEQHAXWRN "${SEQHAXWRN} -Wmissing-prototypes -Wwrite-strings -Wredundant-decls")
SET(SEQHAXWRN "${SEQHAXWRN} -Wchar-subscripts -Wcomment -Wformat=2 -Wwrite-strings")
SET(SEQHAXWRN "${SEQHAXWRN} -Wmissing-declarations -Wredundant-decls -Wnested-externs")
SET(SEQHAXWRN "${SEQHAXWRN} -Wbad-function-cast -Wswitch-enum -Winit-self")
SET(SEQHAXWRN "${SEQHAXWRN} -Wmissing-field-initializers -Wdeclaration-after-statement")
SET(SEQHAXWRN "${SEQHAXWRN} -Wold-style-definition -Waddress -Wmissing-noreturn ")
SET(SEQHAXWRN "${SEQHAXWRN} -Wstrict-overflow=1 -Wextra -Warray-bounds -Wall")

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99 ${SEQHAX_DEPENDS_CFLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SEQHAXWRN}")
SET(CMAKE_C_FLAGS_DEBUG "-ggdb -O0")
SET(CMAKE_C_FLAGS_RELEASE "-O3 -march=native")
include_directories(${CMAKE_SOURCE_DIR}/src
	${CMAKE_BINARY_DIR}
	${CMAKE_SOURCE_DIR}/libqes/src)
link_directories(${CMAKE_BINARY_DIR}/lib)

add_subdirectory(libqes)
add_subdirectory(src)
