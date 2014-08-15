# Targets
ADD_LIBRARY(libseqhax STATIC libseqhax.c ssw.c)
TARGET_LINK_LIBRARIES(libseqhax ${SEQHAX_LIB_DEPENDS})
SET_TARGET_PROPERTIES(libseqhax PROPERTIES OUTPUT_NAME seqhax)

# SSWord - adapter remover
ADD_EXECUTABLE(ssword ssword.c)
TARGET_LINK_LIBRARIES(ssword ${SEQHAX_LIB_DEPENDS} libseqhax)
INSTALL(TARGETS ssword DESTINATION "bin")

