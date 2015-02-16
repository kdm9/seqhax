# Targets
ADD_LIBRARY(libseqhax STATIC seqhax.c ssw.c)
TARGET_LINK_LIBRARIES(libseqhax qes ${SEQHAX_LIB_DEPENDS})
SET_TARGET_PROPERTIES(libseqhax PROPERTIES OUTPUT_NAME seqhax)

# fastq_remove_prefix - Remove a unique 16-mer to each different sample from a lane
ADD_EXECUTABLE(fastq_remove_prefix fastq_remove_prefix.c)
TARGET_LINK_LIBRARIES(fastq_remove_prefix libseqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS fastq_remove_prefix DESTINATION "bin")

# fastq_add_prefix - Add a unique 16-mer to each different sample from a lane
ADD_EXECUTABLE(fastq_add_prefix fastq_add_prefix.c)
TARGET_LINK_LIBRARIES(fastq_add_prefix libseqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS fastq_add_prefix DESTINATION "bin")

# fastq_check - check a fastq file
ADD_EXECUTABLE(fastq_check fastq_check.c)
TARGET_LINK_LIBRARIES(fastq_check libseqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS fastq_check DESTINATION "bin")

# trimadap.c - Heng Li's trimadap, with some fixes he won't merge
ADD_EXECUTABLE(trimadap trimadap.c ksw.c)
TARGET_LINK_LIBRARIES(trimadap ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS trimadap DESTINATION "bin")
