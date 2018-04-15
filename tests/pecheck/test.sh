if [ "${VERBOSE:-no}" == yes ]; then set -x; fi
echo "Testing PECHECK"
rm -rf o && mkdir o

#############
#  Prepare  #
#############
cp il.fq o/expect.fq

#############
#  r1 + r2  #
#############
seqhax pecheck r1.fq r2.fq -o o/tmp.fq >o/stdout 2>o/stderr|| exit 1
cmp o/tmp.fq o/expect.fq || exit 1
grep -q OK o/stdout || exit 1

########
#  il  #
########
seqhax pecheck -i il.fq -o o/tmp.fq  >o/stdout 2>o/stderr|| exit 1
cmp o/tmp.fq o/expect.fq || exit 1
grep -q OK o/stdout || exit 1

#############
#  bad one  #
#############
seqhax pecheck r1.fq r1.fq  >o/stdout 2>o/stderr && exit 1

# Done!
echo "ALL PASSED"
