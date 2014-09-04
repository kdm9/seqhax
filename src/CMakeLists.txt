# Targets
#ADD_LIBRARY(libseqhax STATIC libseqhax.c ssw.c)
#TARGET_LINK_LIBRARIES(libseqhax ${SEQHAX_LIB_DEPENDS})
#SET_TARGET_PROPERTIES(libseqhax PROPERTIES OUTPUT_NAME seqhax)


# fastq_check - check a fastq file
ADD_EXECUTABLE(fastq_check fastq_check.c)
TARGET_LINK_LIBRARIES(fastq_check ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS fastq_check DESTINATION "bin")
