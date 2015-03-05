# Targets
ADD_LIBRARY(libseqhax STATIC seqhax.c ssw.c)
TARGET_LINK_LIBRARIES(libseqhax qes_static ${SEQHAX_LIB_DEPENDS})
SET_TARGET_PROPERTIES(libseqhax PROPERTIES OUTPUT_NAME seqhax)

# fastq_add_prefix - Add a unique 16-mer to each different sample from a lane
ADD_EXECUTABLE(fastq-add-prefix fastq-add-prefix.c)
TARGET_LINK_LIBRARIES(fastq-add-prefix libseqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS fastq-add-prefix DESTINATION "bin")

# fastq_remove_prefix - Remove a unique 16-mer to each different sample from a lane
ADD_EXECUTABLE(fastq-remove-prefix fastq-remove-prefix.c)
TARGET_LINK_LIBRARIES(fastq-remove-prefix libseqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS fastq-remove-prefix DESTINATION "bin")

# fastq_check - check a fastq file
ADD_EXECUTABLE(fastq-check fastq-check.c)
TARGET_LINK_LIBRARIES(fastq-check libseqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS fastq-check DESTINATION "bin")

# ilfastq_len_filter - removes short pairs from an interleaved fastq
ADD_EXECUTABLE(ilfastq-len-filter ilfastq-len-filter.c)
TARGET_LINK_LIBRARIES(ilfastq-len-filter libseqhax ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS ilfastq-len-filter DESTINATION "bin")

# trimadap.c - Heng Li's trimadap, with some fixes he won't merge
ADD_EXECUTABLE(trimadap trimadap.c ksw.c)
TARGET_LINK_LIBRARIES(trimadap ${SEQHAX_LIB_DEPENDS})
INSTALL(TARGETS trimadap DESTINATION "bin")
