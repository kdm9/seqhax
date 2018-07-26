SET(LIBQES_AS_SUBMODULE True)
ADD_SUBDIRECTORY(ext/libqes)

INCLUDE_DIRECTORIES(ext)
INCLUDE_DIRECTORIES(ext/gzstream)
SET(EXT_SOURCES  ext/gzstream/gzstream.C)

file(GLOB CMD_SOURCES progs/*.c progs/*.cc)


# Targets
ADD_EXECUTABLE(seqhax seqhax-main.cc ${CMD_SOURCES} ${EXT_SOURCES})
TARGET_LINK_LIBRARIES(seqhax qes_static ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS seqhax DESTINATION "bin")
