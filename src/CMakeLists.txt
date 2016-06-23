file(GLOB LIBQES_SOURCES libqes/qes_*.c)
# Targets
ADD_EXECUTABLE(seqhax seqhax.c seqhax-main.c ${LIBQES_SOURCES})
TARGET_LINK_LIBRARIES(seqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS seqhax DESTINATION "bin")
