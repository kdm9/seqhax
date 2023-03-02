SET(LIBQES_AS_SUBMODULE True)
ADD_SUBDIRECTORY(lib/libqes)

INCLUDE_DIRECTORIES(lib)
INCLUDE_DIRECTORIES(lib/gzstream)
SET(EXT_SOURCES  lib/gzstream/gzstream.C)


# Targets
file(GLOB CMD_SOURCES seqhax/progs/*.c seqhax/progs/*.cc)
ADD_EXECUTABLE(seqhax seqhax/seqhax-main.cc ${CMD_SOURCES} ${EXT_SOURCES})
TARGET_LINK_LIBRARIES(seqhax qes_static ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS seqhax DESTINATION "bin")

ADD_EXECUTABLE(htshax htshax/main.cc htshax/htshax_bcfhist.cc ${EXT_SOURCES})
TARGET_LINK_LIBRARIES(htshax qes_static ${SEQHAX_LIB_DEPENDS} ${HTSLIB_LIBRARIES})
INSTALL(TARGETS htshax DESTINATION "bin")
