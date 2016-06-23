file(GLOB LIB_SOURCES libqes/qes_*.c lib/*.c)
file(GLOB CMD_SOURCES progs/*.c)

include_directories(libs)
# Targets
ADD_EXECUTABLE(seqhax seqhax-main.c ${LIB_SOURCES} ${CMD_SOURCES})
TARGET_LINK_LIBRARIES(seqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS seqhax DESTINATION "bin")
